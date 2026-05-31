# Container format — deep dive

The `.pack` is PRST's optional bundle format: instead of scattering a proof's many point and product files across the directory, `-proof … pack` writes them all into one append-friendly container. It's the *implementation* behind a contract two other docs only sketch — `proof-system.md` §6 (the `.pack` option) and `state-serialization.md` §8 (`FilePacked` as the on-disk path). This doc covers how that container is actually laid out and read.

**It is a live path** (grep-confirmed, not dormant): `-proof … pack` constructs a `container::FileContainer` (`boinc.cpp:315`, `prst.cpp:366`), and the proof system writes its points/products as `FilePacked` files into it (`proof.cpp:189,199`; `FilePacked`/`FilePackedWriter` in `file.cpp:385-509`). It's *optional* — a normal `-proof save/build` uses individual `File`s — but when `pack` is requested, this code runs.

> **Scope.** `container.cpp` is ~1,897 lines, including a hand-rolled JSON parser, a base64 codec, and corruption recovery. The class model (from `container.h`) and the **on-disk layout** (from the verified `FileContainer::open` parser) are documented here; the byte-level *writer* framing, the codec internals, and the exact recovery walk are pointed to, not reproduced. External-format risk handled the same way as `curves-and-polynomials.md`: claims are grounded in the header + the parser I read, not inferred from the bodies I didn't.

Source files: `framework/container.h`, `framework/container.cpp`. The `File`-layer bridge: `framework/file.cpp:385-509` (`FilePacked`, `FilePackedWriter`).

Prereqs / companions: `state-serialization.md` (`FilePacked` is a `File` subclass; this is its backend, and a packed file's *payload* is still the `Writer`/`Reader` framing from that doc), `proof-system.md` §6 (the `-proof pack` pipeline — the only consumer).

## 1. The on-disk format

A `.pack` is a sequence of **newline-delimited JSON-array records** (`FileContainer::open`, `container.cpp:1681-1796`):

```
["PK", 1, …]                         ← header line (≤32 chars): magic + format version 1
[<chunk_size>, <stream_id>, <codec>] ← chunk record header (≤1023 chars), then <chunk_size> raw bytes
[<chunk_size>, <stream_id>, <codec>]
…
[0, 0]                               ← terminator (a stream-0 chunk of size 0)
```

- **Header.** `root[0]` is a magic string — `"PK"` (zip-like), `"7f"`, or `"\xD0\x9A"` — and `root[1]` must be `1` (the format version). Anything else ⇒ `EMPTY`.
- **Chunk records.** Each is a JSON line `[chunk_size, stream_id, optional-codec]` followed by exactly `chunk_size` bytes:
  - **`stream_id == 0` → an index chunk.** The `chunk_size` bytes are themselves a JSON array of file entries `{file, stream, size, offset, md5}`, which populate the `Index`. A size-0 stream-0 chunk is the **terminator**.
  - **`stream_id > 0` → a data chunk** for that logical stream. The parser records the chunk's `{offset, size, codec}` and *skips* the bytes (lazy — data is read only when a file is requested). An optional `codec` JSON at `root[2]` says how to decode the chunk (e.g. base64).
- **Files span streams.** A file's bytes are the concatenation of its stream's data chunks, sliced to `[offset, offset+size)`; `md5` (when present) verifies it. One stream can hold several files (sorted by offset).

Records are **CRLF-terminated** (`\r\n`). The writer side is symmetric and verified: `ChunkedWriteStream::write_chunk` (`container.cpp:834-852`) emits exactly `[buffer_size, stream_id, codec?]\r\n` followed by the raw bytes — the same record `open` parses. Data is buffered up to `MAX_BUFFER_SIZE` (~16 MB) per chunk; the index chunk and the `[0,0]` terminator are emitted by the `Packer` on close (the per-stream terminator in `ChunkedWriteStream::close` is commented out, `:864-865`).

So the format interleaves data chunks with periodic index chunks. **It's append-friendly**: to add files you append more data chunks and a fresh index chunk at the end; on open, a later index entry for a filename **overrides** the earlier one (`open` does `entry.reset()` then reassigns, `container.cpp:1743-1746`). That's how a proof accumulates points into one `.pack` incrementally.

## 2. `FileContainer::open` — the central method

The parser is the authoritative definition of the format. Annotated (`container.cpp:1681-1796`):

```cpp
int FileContainer::open()
{
    // header: a JSON array ["PK"|"7f"|"\xD0\x9A", 1, …]; else EMPTY
    if (!_stream->readline(st, 32) || !json.parse(st) || !json.root().is_array()
        || json.root()[0] not a known magic || json.root()[1] != 1) return EMPTY;

    while (true) {
        // each chunk: [chunk_size, stream_id, codec?]
        if (!_stream->readline(st, 1023) || !json … size<2 …) return UNEXPECTED_END;
        int64_t chunk_size = json.root()[0].value_int();
        int64_t next = _stream->position() + chunk_size;
        if (next > _stream->length()) return UNEXPECTED_END;
        int64_t stream_id = json.root()[1].value_int();

        if (stream_id == 0) {                          // index chunk
            if (chunk_size == 0) break;                //   size 0 = terminator
            read chunk_size bytes; json.parse(them);   //   a JSON array of file entries
            for (file : entries) {                     //   {file, stream, size, offset, md5}
                auto& entry = _index[file["file"]];     //   later entry OVERRIDES earlier
                entry.stream = file["stream"]; entry.size/offset/md5 = …;
            }
        } else {                                       // data chunk for stream_id
            auto& stream = _streams[stream_id];
            stream.chunks.push_back({offset: pos, size: chunk_size, codec: root[2]});
            _stream->set_position(next);               //   skip the bytes (lazy)
        }
    }

    // validate: per stream, sort files by offset; flag inconsistent ones → on_corrupted, INDEX_CORRUPTED
    // offset > 16777215 (one ~16 MB chunk) ⇒ INDEX_CORRUPTED
    if (_stream->position() < _stream->length()) return EXCESS_DATA;
    return error;
}
```

The return is a `container_error` code: `OK` / `EMPTY` (no/blank/bad-header) / `UNEXPECTED_END` (truncated) / `INDEX_CORRUPTED` (inconsistent entries — `on_corrupted` moves them to `_corrupted` and parsing continues) / `EXCESS_DATA` (trailing bytes past the last chunk). The caller treats `OK` and `EMPTY` as usable (an empty container is fine to append to — `boinc.cpp:316`), the rest as corruption to warn about.

## 3. Class model & reference

Three layers in `container.h`:

**JSON** (`:23-116`) — a from-scratch JSON parser/serializer. A `Node` is a `variant<null, value-buffer, array, object>`; objects are a *sorted* `(keys, values)` flat-map (binary-searched `operator[]`/`contains`). Used for every record line and the index. (It's here so the container has no external JSON dependency.)

**Streams** (`:120-339`) — a `ReadStream`/`WriteStream` abstraction with concretes `MemoryStream`, `LimitedReadStream`, `FileStream`, plus the format-specific decorators:
- `ChunkedWriteStream` — buffers up to `MAX_BUFFER_SIZE` (16,777,200 ≈ 16 MB) then emits a data chunk (`write_chunk`); carries a per-chunk `codec_json`.
- `Base64CoderStream` / `Base64DecoderReadStream` — the base64 codec a chunk's `codec` field selects.

**Container** (`:222-489`):

| Type | Role |
|---|---|
| `FileDesc` | one index entry: `{name, size, stream, offset, md5, data}` |
| `Index` | a `set<FileDesc>` sorted by name (`contains`/`operator[]` by filename) |
| `Packer` (+ `Packer::Writer`) | **writes** a container: `add_writer()` → `Writer`; `Writer::add_file(name, …)` → a `WriteStream`; `add_codec`; `md5` flag. Built from a filename, a `WriteStream`, or a `FileContainer`. |
| `FileContainer` (+ inner `Reader`) | **reads/appends** a container: `open()` parses it; `read_file(filename)` → a `FileDesc` whose `data` `ReadStream` reassembles the stream's chunks and applies codecs; `reopen`/`close`. Holds `_index`, `_streams` (id → `ChunkStream{chunks, files}`), `_corrupted`. |
| `Unpacker` | a `WriteStream` that dispatches extracted files via an `on_file` callback (the inverse of `Packer`). |

The `File`-layer bridge (`file.cpp`): `FilePacked` is a `File` whose `read_buffer` pulls bytes from `FileContainer::read_file`, and whose `get_writer` returns a `FilePackedWriter` wrapping a `Packer` (`file.cpp:412-487`). So the rest of PRST writes a packed proof point exactly like any other `File` — the container is hidden behind the `File` interface (`state-serialization.md` §3).

## 4. Lifecycle

**Write** (`-proof … pack`): `Proof::init_files` wraps each point/product file as a `FilePacked` over the shared `FileContainer` (`proof.cpp:189,199`). Each `File::write` opens a `FilePackedWriter` → a `Packer::Writer` → an `add_file` `WriteStream` (a `ChunkedWriteStream`, optionally base64-coded), writes the framed payload (the same magic/`TYPE`/fingerprint record from `state-serialization.md`), and on close flushes data chunk(s) + updates the index. The container is `close()`d at the end of the run (`fermat.cpp:397-398`).

**Read / resume**: constructing the `FileContainer` runs `open()` (§2), building the in-memory `Index` and `_streams`. `read_file(name)` then hands back a `Reader` that lazily concatenates that file's chunks (decoding per the chunk codec) into the `FilePacked` buffer, which the normal `Reader` framing then parses.

**Append**: re-opening for write and adding more files appends new data + index chunks; the override rule (§1) means the newest index wins. This is what lets a long proof grow its `.pack` over many checkpoints.

## 5. Pitfalls

- **`EMPTY` is not an error.** An absent or zero-length or blank-header file returns `container_error::EMPTY`, and callers treat that as "fresh container, fine to append" (`boinc.cpp:316` warns only on errors *other* than `OK`/`EMPTY`). Don't treat `EMPTY` as failure.
- **Later index entries override earlier ones.** A filename appearing in two index chunks resolves to the *last* (`open` resets and reassigns). That's the append mechanism — but it also means a truncated/duplicated append can silently shadow good data; the offset/size validation and MD5 are the guards.
- **Offsets are bounded by the 16 MB chunk size.** `open` flags `offset > 16777215` as `INDEX_CORRUPTED` — a file's offset is *within* its stream's chunk addressing, not a global file position. Reading the offset as an absolute byte position into the `.pack` is wrong.
- **Data is read lazily and skipped on open.** `open` records chunk positions and seeks past the bytes; the payload is only materialized on `read_file`. A consumer that assumes `open` validated every byte (beyond the index) is mistaken — MD5 is checked at read time, not open time.
- **The payload inside a packed file is still the `state-serialization.md` framing.** The container is a *bundling* layer; each entry's bytes are an ordinary magic+appid+`TYPE`+fingerprint record. A `.pack` bug and a checkpoint-format bug are different layers — check which one you're in.
- **It's a bespoke format with a hand-rolled JSON parser.** Not zip (despite the `"PK"` magic) and not a standard archive — don't point external tools at it. The parser, base64 codec, and recovery are all in `container.cpp`.

## 6. Quick reference

On-disk: newline-delimited JSON records. Header `["PK"|"7f"|"\xD0\x9A", 1, …]`; then `[size, stream_id, codec?]` + `size` bytes per chunk; `stream_id 0` = index (JSON array of `{file,stream,size,offset,md5}`), `>0` = data; `[0,0]` terminates.

| You want… | Where |
|---|---|
| Read a packed file | `FileContainer::read_file(name)` → `FileDesc::data` |
| Write into a container | `Packer::add_writer()` → `Writer::add_file(...)` |
| The `File`-API bridge | `FilePacked` / `FilePackedWriter` (`file.cpp`) |
| Enable the `.pack` | `-proof … pack [name]` (`prst.cpp:366`, `boinc.cpp:315`) |
| Interpret an open error | `container_error`: `OK`/`EMPTY` usable; `UNEXPECTED_END`/`INDEX_CORRUPTED`/`EXCESS_DATA` = trouble |
| Add a codec to a chunk | `Writer::add_codec` / `ChunkedWriteStream` `codec_json` (base64) |

## 7. Open questions / non-coverage

- **The `Packer`-level write path.** The *data*-chunk framing is verified (`open` for the reader, `write_chunk` `:834-852` for the writer — they're symmetric). What's *not* traced: how `Packer`/`Packer::Writer::close` emit the **index chunk** and the `[0,0]` terminator, the `add_file` → stream wiring, and flush ordering across writers. Read `Packer` before changing how the index is written.
- **Codec internals.** The base64 coder/decoder (`Base64CoderStream`/`Base64DecoderReadStream`) and the `codec` JSON schema beyond "base64 selects it" aren't dissected; whether other codecs exist or are planned is unexamined.
- **Corruption-recovery behavior.** `on_corrupted` (`container.cpp:1798`) moves bad entries to `_corrupted` and parsing continues, but the full recovery walk (what's salvageable, how `_cur_reader`/`_cur_file` interact) was not traced.
- **The JSON parser.** `JSON::Node::parse` is a hand-rolled parser; its exact grammar/limits (number formats, escaping, the `value-buffer` lazy-parse) are out of scope — treated as "parses the records," not audited.
- **Concurrency / partial-write durability.** Unlike the local `File` write (atomic `.new`+rename, `state-serialization.md` §5), the container is appended in place; what guarantees hold if the process dies mid-append (beyond the `open`-time corruption codes) wasn't analyzed.
