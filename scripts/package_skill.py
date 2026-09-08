"""Package only the explicitly reviewed files in release-files.txt."""
from pathlib import Path
import hashlib
import zipfile
import yaml


def release_files(root):
    root = Path(root).resolve()
    files = []
    for name in (root / "release-files.txt").read_text(encoding="utf-8").splitlines():
        if not name or name.startswith("#"):
            continue
        relative = Path(name)
        if relative.is_absolute() or relative.drive or ".." in relative.parts:
            raise SystemExit(f"Invalid release path: {name}")
        path = root / relative
        if path.resolve() != path or not path.is_file():
            raise SystemExit(f"Missing or linked release file: {name}")
        if path in files:
            raise SystemExit(f"Duplicate release path: {name}")
        files.append(path)
    if not files:
        raise SystemExit("Empty release manifest")
    return files


def package(root, output_dir=None):
    root = Path(root).resolve()
    files = release_files(root)
    metadata = yaml.safe_load((root / "SKILL.md").read_text(encoding="utf-8").split("---", 2)[1])
    version = metadata["metadata"]["version"]
    output = Path(output_dir or root.parent) / f"{root.name}-{version}.zip"
    sidecar = output.with_suffix(".sha256.txt")
    if output.exists() or sidecar.exists():
        raise SystemExit("Versioned package already exists")
    with zipfile.ZipFile(output, "x", zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, (Path(root.name) / path.relative_to(root)).as_posix())
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise SystemExit("Release archive integrity check failed")
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    sidecar.write_text(f"name={root.name}\nversion={version}\nzip={output.name}\nn_files={len(files)}\nbytes={output.stat().st_size}\nsha256={digest}\n", encoding="utf-8")
    print(output)
    return output


if __name__ == "__main__":
    package(Path(__file__).resolve().parents[1])

