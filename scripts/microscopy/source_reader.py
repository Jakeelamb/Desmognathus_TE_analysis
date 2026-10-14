"""Optional pinned Bio-Formats metadata and bounded, unchanged green-channel crops."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from xml.etree import ElementTree

# Import the existing environment guard and hashing; do not alter its toolchain.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from recovery_backend import active_prefix, clean_native_env, no_symlink, sha256

VERSION = "8.5.0"
JAR_SHA256 = "c6e60665d53a334b66e4d635340151f403dfe57a64704c573dd4c03b873befb9"
JAR_URL = (
    f"https://downloads.openmicroscopy.org/bio-formats/{VERSION}/artifacts/bioformats_package.jar"
)
ROOT = Path(__file__).resolve().parents[2]
MAX_PIXELS = 1_000_000
MAX_TILE_PIXELS = 4 * MAX_PIXELS
MAX_JAR_BYTES = 150_000_000
TIMEOUT = 300


def now() -> str:
    return datetime.now(UTC).isoformat()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def cache_path(root: Path) -> Path:
    return root / ".pixi/source-reader-tools" / f"{VERSION}-{JAR_SHA256}"


def java_path(root: Path) -> Path:
    prefix = active_prefix(root)
    java = (prefix / "bin/java").resolve()
    if not java.is_file() or not java.is_relative_to(prefix):
        raise ValueError("Java must belong to the locked recovery environment")
    return java


def tool_pins(root: Path, java: Path) -> dict:
    home = java.parent.parent
    paths = [java, home / "release", home / "lib/modules", home / "lib/server/libjvm.so"]
    return {
        "bioformats_version": VERSION,
        "jar_url": JAR_URL,
        "jar_sha256": JAR_SHA256,
        "pixi_lock_sha256": sha256(root / "pixi.lock"),
        "reader_sha256": sha256(Path(__file__)),
        "java_pins": {str(path): sha256(path) for path in paths},
    }


def setup(root: Path, seed: Path | None = None) -> Path:
    java = java_path(root)
    cache = cache_path(root)
    no_symlink(cache)
    if cache.exists():
        return load_jar(root)
    cache.mkdir(parents=True, exist_ok=False)
    receipt = {
        "start_utc": now(),
        "status": "failed",
        "jar_url": JAR_URL,
        "jar_sha256_expected": JAR_SHA256,
        "seed": str(seed) if seed else None,
    }
    try:
        receipt["tool_pins"] = tool_pins(root, java)
        destination = cache / "bioformats_package.jar"
        if seed is not None:
            with seed.open("rb") as source:
                _copy_bounded(source, destination)
        else:
            with urllib.request.urlopen(JAR_URL, timeout=60) as source:
                _copy_bounded(source, destination)
        receipt["jar_sha256_actual"] = sha256(destination)
        if receipt["jar_sha256_actual"] != JAR_SHA256:
            raise ValueError("Bio-Formats JAR checksum mismatch")
        receipt["status"] = "complete"
    except Exception as error:
        receipt["error"] = repr(error)
        raise
    finally:
        receipt["end_utc"] = now()
        write_json(cache / "setup.json", receipt)
    return destination


def _copy_bounded(source, destination: Path) -> None:
    count = 0
    with destination.open("xb") as target:
        while block := source.read(1024 * 1024):
            count += len(block)
            if count > MAX_JAR_BYTES:
                raise ValueError("Bio-Formats download exceeds size bound")
            target.write(block)


def load_jar(root: Path) -> Path:
    cache = cache_path(root)
    no_symlink(cache)
    receipt = cache / "setup.json"
    jar = cache / "bioformats_package.jar"
    no_symlink(receipt)
    no_symlink(jar)
    if not receipt.is_file() or not jar.is_file():
        raise ValueError("Optional source reader missing; run explicit source-reader setup")
    if json.loads(receipt.read_text()).get("status") != "complete":
        raise ValueError("Failed source reader cache retained; no in-place repair permitted")
    if sha256(jar) != JAR_SHA256:
        raise ValueError("Source reader cache changed; no in-place repair permitted")
    return jar


def crop_bounds(x: int, y: int, width: int, height: int) -> None:
    if any(type(v) is not int for v in (x, y, width, height)):
        raise ValueError("Crop coordinates and dimensions must be integers")
    if x < 0 or y < 0 or width < 1 or height < 1 or width * height > MAX_PIXELS:
        raise ValueError(
            "Crop requires nonnegative coordinates, positive dimensions and <=1M pixels"
        )


def fresh_run(root: Path, out: Path, source: Path) -> Path:
    out = out.absolute()
    no_symlink(out)
    runs = (root / ".pixi/runs").resolve()
    resolved = out.resolve()
    if not resolved.is_relative_to(runs) or resolved == runs:
        raise ValueError("Outputs must be fresh directories below .pixi/runs/")
    if resolved.is_relative_to(source.resolve().parent):
        raise ValueError("Output in source or companion tree refused")
    out.mkdir(parents=True, exist_ok=False)
    return out


def source_pin(source: Path, hash_mode: str = "full") -> dict:
    if hash_mode not in {"full", "stat-only"}:
        raise ValueError("Source hash mode must be full or stat-only")
    if not source.is_file() or source.suffix.lower() not in {".vsi", ".tif", ".tiff", ".btf"}:
        raise ValueError("Source must be an existing VSI or TIFF file")
    stat = source.stat()
    return {
        "path": str(source.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "sha256": sha256(source) if hash_mode == "full" else None,
        "full_payload_hash_computed": hash_mode == "full",
        "hash_scope": (
            "named source file only; companion ETS files not hashed"
            if hash_mode == "full"
            else "filesystem stat only; source payload not hashed"
        ),
        "access_policy": "read operations only; source mount is not declared read-only",
    }


def run_java(argv: list[str], out: Path) -> dict:
    command = {"argv": argv, "start_utc": now(), "exit_code": None}
    stdout, stderr = out / "command.stdout.log", out / "command.stderr.log"
    try:
        with stdout.open("xb") as output, stderr.open("xb") as errors:
            result = subprocess.run(
                argv,
                cwd=out,
                env=clean_native_env(),
                stdin=subprocess.DEVNULL,
                stdout=output,
                stderr=errors,
                timeout=TIMEOUT,
                check=False,
            )
        command["exit_code"] = result.returncode
    except (OSError, subprocess.TimeoutExpired) as error:
        command["error"] = repr(error)
    command.update(
        end_utc=now(),
        stdout_sha256=sha256(stdout),
        stderr_sha256=sha256(stderr),
        stdout_path=str(stdout),
        stderr_path=str(stderr),
    )
    write_json(out / "command.json", command)
    return command


def validate_crop(path: Path, width: int, height: int):
    import numpy as np
    import tifffile

    with tifffile.TiffFile(path) as image:
        if (
            len(image.pages) != 1
            or len(image.series) != 1
            or len(image.series[0].levels) != 1
            or image.pages[0].subifds
        ):
            raise ValueError("Crop must contain exactly one TIFF plane")
        page = image.pages[0]
        if page.samplesperpixel != 1 or page.photometric.name != "MINISBLACK":
            raise ValueError("Crop must be a single grayscale MINISBLACK plane")
        if page.dtype != np.dtype("uint8") or page.shape != (height, width):
            raise ValueError("Crop must be uint8 with exact requested dimensions")
        pixels = page.asarray()
        if pixels.dtype != np.uint8 or pixels.shape != (height, width):
            raise ValueError("Decoded crop pixels have invalid dtype or shape")
    return pixels


def tiled_patch(page, source: Path, x: int, y: int, width: int, height: int, operation: dict):
    """Decode only intersecting tiles through tifffile's canonical segment decoder."""
    import numpy as np

    tile_width, tile_height = page.tilewidth, page.tilelength
    columns = (page.imagewidth + tile_width - 1) // tile_width
    rows = (page.imagelength + tile_height - 1) // tile_height
    expected_bytes = tile_width * tile_height * 3
    first_col, last_col = x // tile_width, (x + width - 1) // tile_width
    first_row, last_row = y // tile_height, (y + height - 1) // tile_height
    count = (last_col - first_col + 1) * (last_row - first_row + 1)
    operation.update(planned_tile_count=count, decoded_tile_count=0, encoded_bytes_read=0)
    if count * tile_width * tile_height > MAX_TILE_PIXELS:
        raise ValueError("Intersecting tiles exceed the 4M decoded-pixel read bound")
    if len(page.dataoffsets) != columns * rows or len(page.databytecounts) != columns * rows:
        raise ValueError("Tile directory does not match actual TIFF page geometry")
    patch = np.empty((height, width), dtype=np.uint8)
    reads = []
    operation["tile_reads"] = reads
    size = source.stat().st_size
    with source.open("rb") as handle:
        for row in range(first_row, last_row + 1):
            for col in range(first_col, last_col + 1):
                index = row * columns + col
                offset, byte_count = page.dataoffsets[index], page.databytecounts[index]
                if offset < 1 or byte_count != expected_bytes or offset + byte_count > size:
                    raise ValueError(
                        "Tile byte range is missing, truncated or inconsistent with RGB8"
                    )
                handle.seek(offset)
                encoded = handle.read(byte_count)
                reads.append(
                    {
                        "tile_index": index,
                        "offset": offset,
                        "byte_count": byte_count,
                        "actual_bytes_read": len(encoded),
                    }
                )
                operation["encoded_bytes_read"] += len(encoded)
                if len(encoded) != byte_count:
                    raise ValueError("Truncated TIFF tile read")
                decoded, position, shape = page.decode(encoded, index)
                top, left = row * tile_height, col * tile_width
                if (
                    decoded is None
                    or decoded.dtype != np.uint8
                    or decoded.shape != (1, tile_height, tile_width, 3)
                    or tuple(shape) != decoded.shape
                    or tuple(position) != (0, 0, top, left, 0)
                ):
                    raise ValueError(
                        "Canonical TIFF tile decoder returned unexpected layout or dtype"
                    )
                y0, y1 = max(y, top), min(y + height, top + tile_height)
                x0, x1 = max(x, left), min(x + width, left + tile_width)
                patch[y0 - y : y1 - y, x0 - x : x1 - x] = decoded[
                    0, y0 - top : y1 - top, x0 - left : x1 - left, 1
                ]
                operation["decoded_tile_count"] += 1
    return patch


def tiff_page_read(
    source: Path, out: Path, action: str, x: int, y: int, width: int, height: int
) -> tuple[Path, dict]:
    """Read actual page 0 headers and bounded pixels; never interpret OME series."""
    import numpy as np
    import tifffile

    operation = {
        "execution_kind": "in-process Python; no native subprocess",
        "method": "tifffile actual page-0 header",
        "start_utc": now(),
        "status": "failed",
        "process_exit_code": None,
        "pixel_scope": "TIFF page 0; embedded OME series ignored",
    }
    stdout, stderr = out / "command.stdout.log", out / "command.stderr.log"
    stdout.write_text("")
    stderr.write_text("")
    try:
        with tifffile.TiffFile(source, mode="r") as image:
            if len(image.pages) != 1:
                raise ValueError("TIFF-page reader requires exactly one actual TIFF page")
            page = image.pages[0]
            header = {
                "page": 0,
                "page_count": len(image.pages),
                "bigtiff": image.is_bigtiff,
                "width": page.imagewidth,
                "height": page.imagelength,
                "shape": list(page.shape),
                "dtype": str(page.dtype),
                "samples_per_pixel": page.samplesperpixel,
                "bits_per_sample": page.bitspersample,
                "photometric": page.photometric.name,
                "planar_config": page.planarconfig.name,
                "compression": page.compression.name,
                "extra_samples": list(page.extrasamples),
                "subifds": page.subifds,
                "memmappable": page.is_memmappable,
                "tiled": page.is_tiled,
                "tile_width": page.tilewidth,
                "tile_height": page.tilelength,
                "orientation": int(page.tags["Orientation"].value)
                if "Orientation" in page.tags
                else 1,
                "predictor": int(page.predictor),
                "fill_order": int(page.fillorder),
                "data_offsets": list(page.dataoffsets),
                "data_byte_counts": list(page.databytecounts),
                "coordinate_authority": "actual TIFF page header",
                "embedded_ome_series_validity": "not evaluated; ignored for pixel indexing",
            }
            operation["header"] = header
            stdout.write_text(json.dumps(header, indent=2) + "\n")
            description = out / "embedded-description.txt"
            description.write_text(page.description or "")
            header["embedded_description"] = {
                "path": str(description),
                "sha256": sha256(description),
                "contains_ome_text": "<OME" in (page.description or ""),
                "scope": "embedded ImageDescription text, retained without OME-series validation",
            }
            if (
                page.subifds
                or page.compression != 1
                or (not page.is_tiled and not page.is_memmappable)
                or page.photometric.name != "RGB"
                or page.planarconfig.name != "CONTIG"
                or page.samplesperpixel != 3
                or page.extrasamples
                or page.dtype != np.dtype("uint8")
                or page.bitspersample != 8
                or page.shape != (page.imagelength, page.imagewidth, 3)
                or header["orientation"] != 1
                or page.predictor != 1
                or page.fillorder != 1
                or (
                    page.is_tiled
                    and (page.tiledepth != 1 or page.tilewidth * page.tilelength > MAX_PIXELS)
                )
            ):
                raise ValueError(
                    "TIFF-page reader requires uncompressed tiled or memmappable CONTIG RGB3 uint8 "
                    "with orientation 1, no predictor, extra samples or subresolutions"
                )
            if action == "crop" and (x + width > page.imagewidth or y + height > page.imagelength):
                raise ValueError("Crop exceeds actual TIFF page header bounds")
            patch = None
            if action == "crop" and page.is_tiled:
                operation["method"] = (
                    "tifffile canonical decode of intersecting page-0 tile byte ranges"
                )
                patch = tiled_patch(page, source, x, y, width, height, operation)
        if action == "metadata":
            asset = out / "metadata.json"
            write_json(asset, header)
        else:
            if patch is None:
                operation["method"] = "tifffile read-only page-0 memmap slice"
                mapped = tifffile.memmap(source, page=0, mode="r")
                try:
                    if mapped.dtype != np.uint8 or mapped.shape != tuple(header["shape"]):
                        raise ValueError(
                            "Memory-map shape or dtype differs from checked page header"
                        )
                    if mapped.flags.writeable:
                        raise ValueError("Source memory-map must be read-only")
                    patch = mapped[y : y + height, x : x + width, 1].copy()
                finally:
                    mapped._mmap.close()
            asset = out / "crop.tiff"
            tifffile.imwrite(asset, patch, photometric="minisblack", metadata={"axes": "YX"})
        operation["status"] = "complete"
        stdout.write_text(json.dumps(operation, indent=2) + "\n")
    except Exception as error:
        operation["error"] = repr(error)
        stderr.write_text(repr(error) + "\n")
        raise
    finally:
        operation.update(
            end_utc=now(),
            stdout_path=str(stdout),
            stderr_path=str(stderr),
            stdout_sha256=sha256(stdout),
            stderr_sha256=sha256(stderr),
        )
        write_json(out / "command.json", operation)
    return asset, operation


def read_source(
    root: Path,
    action: str,
    source: Path,
    out: Path,
    *,
    x: int = 0,
    y: int = 0,
    width: int = 0,
    height: int = 0,
    reader: str = "bioformats",
    source_hash: str = "full",
) -> Path:
    if action not in {"metadata", "crop"}:
        raise ValueError("Unknown source reader action")
    if reader not in {"bioformats", "tiff-page"}:
        raise ValueError("Unknown source reader method")
    if action == "crop":
        crop_bounds(x, y, width, height)
    source = source.resolve()
    # Validate before writing; pin within the receipt so acquisition failures remain visible.
    if not source.is_file():
        raise ValueError("Source file missing")
    out = fresh_run(root, out, source)
    receipt = {
        "schema": "microscopy-source-reader-v1",
        "action": action,
        "start_utc": now(),
        "status": "failed",
        "source_path": str(source),
        "selection": {"series": 0, "resolution": 0},
        "transforms": [],
        "companion_ets_hashes_computed": False,
        "reader": reader,
        "source_hash_mode": source_hash,
    }
    try:
        receipt["source"] = source_pin(source, source_hash)
        if reader == "tiff-page":
            active_prefix(root)
            import numpy as np
            import tifffile

            receipt["selection"] = {"page": 0, "sample": 1} if action == "crop" else {"page": 0}
            receipt["tools"] = {
                "reader_sha256": sha256(Path(__file__)),
                "pixi_lock_sha256": sha256(root / "pixi.lock"),
                "numpy_version": np.__version__,
                "tifffile_version": tifffile.__version__,
            }
            receipt["coordinate_authority"] = "actual TIFF page header; embedded OME series ignored"
            if action == "crop":
                receipt["crop"] = {"x": x, "y": y, "width": width, "height": height}
            asset, operation = tiff_page_read(source, out, action, x, y, width, height)
            receipt["command"] = operation
            if action == "crop":
                pixels = validate_crop(asset, width, height)
                receipt["validation"] = {
                    "dtype": str(pixels.dtype),
                    "shape": list(pixels.shape),
                    "single_plane_grayscale": True,
                }
            receipt["artifact"] = {"path": str(asset), "sha256": sha256(asset)}
            if source_pin(source, source_hash) != receipt["source"]:
                raise ValueError("Named source changed during read")
            receipt["status"] = "complete"
            return out
        java = java_path(root)
        receipt["tools"] = tool_pins(root, java)
        jar = load_jar(root)
        argv = [
            str(java),
            "-Xmx2g",
            "-Djava.awt.headless=true",
            f"-Djava.io.tmpdir={out}",
            "-cp",
            str(jar),
        ]
        if action == "metadata":
            argv += [
                "loci.formats.tools.ImageInfo",
                "-nopix",
                "-omexml",
                "-no-upgrade",
                "-novalid",
                "-noflat",
                "-series",
                "0",
                str(source),
            ]
        else:
            receipt["selection"].update(channel=1, z=0, t=0)
            receipt["crop"] = {"x": x, "y": y, "width": width, "height": height}
            argv += [
                "loci.formats.tools.ImageConverter",
                "-no-upgrade",
                "-series",
                "0",
                "-separate",
                "-channel",
                "1",
                "-z",
                "0",
                "-timepoint",
                "0",
                "-crop",
                f"{x},{y},{width},{height}",
                "-nooverwrite",
                str(source),
                str(out / "crop.tiff"),
            ]
        command = run_java(argv, out)
        receipt["command"] = command
        if command["exit_code"] != 0:
            raise ValueError("Bio-Formats command failed; inspect retained command logs")
        if action == "metadata":
            raw = (out / "command.stdout.log").read_text()
            start, end = raw.find("<OME"), raw.rfind("</OME>")
            if start < 0 or end < start:
                raise ValueError("Bio-Formats did not emit OME-XML")
            xml = raw[start : end + len("</OME>")]
            tree = ElementTree.fromstring(xml)
            if not tree.findall("{*}Image"):
                raise ValueError("OME-XML contains no image metadata")
            asset = out / "metadata.ome.xml"
            asset.write_text(xml + "\n")
        else:
            asset = out / "crop.tiff"
            pixels = validate_crop(asset, width, height)
            receipt["validation"] = {
                "dtype": str(pixels.dtype),
                "shape": list(pixels.shape),
                "single_plane_grayscale": True,
            }
        receipt["artifact"] = {"path": str(asset), "sha256": sha256(asset)}
        if source_pin(source, source_hash) != receipt["source"]:
            raise ValueError("Named source changed during read")
        receipt["status"] = "complete"
    except Exception as error:
        if "command" not in receipt and (out / "command.json").is_file():
            receipt["command"] = json.loads((out / "command.json").read_text())
        receipt["error"] = repr(error)
        raise
    finally:
        receipt["end_utc"] = now()
        write_json(out / "receipt.json", receipt)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    install = sub.add_parser("setup", help="Explicitly acquire the checksum-pinned optional JAR")
    install.add_argument("--jar", type=Path, help="Use a locally acquired official JAR seed")
    for action in ("metadata", "crop"):
        command = sub.add_parser(action)
        command.add_argument("--source", type=Path, required=True)
        command.add_argument("--out", type=Path, required=True)
        command.add_argument("--reader", choices=("bioformats", "tiff-page"), default="bioformats")
        command.add_argument("--source-hash", choices=("full", "stat-only"), default="full")
        if action == "crop":
            for option in ("x", "y", "width", "height"):
                command.add_argument(f"--{option}", type=int, required=True)
    args = parser.parse_args(argv)
    try:
        if args.action == "setup":
            print(setup(ROOT, args.jar))
        else:
            bounds = (
                {name: getattr(args, name) for name in ("x", "y", "width", "height")}
                if args.action == "crop"
                else {}
            )
            print(
                read_source(
                    ROOT,
                    args.action,
                    args.source,
                    args.out,
                    reader=args.reader,
                    source_hash=args.source_hash,
                    **bounds,
                )
            )
    except (OSError, ValueError, ElementTree.ParseError) as error:
        print(f"source-reader: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
