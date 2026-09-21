#!/usr/bin/env bash
set -euo pipefail

VERSION="${1:-}"
TARGET="${2:-all}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$REPO_ROOT/.release-venv"
VENV_PYTHON="$VENV_DIR/bin/python"

if [[ -z "$VERSION" ]]; then
  echo "Usage: ./release.sh <version> [check|test|pypi|all]"
  echo "Examples:"
  echo "  ./release.sh 0.0.2 check"
  echo "  ./release.sh 0.0.2 test"
  echo "  ./release.sh 0.0.2 pypi"
  echo "  ./release.sh 0.0.2 all"
  exit 1
fi

case "$TARGET" in
  check|test|pypi|all)
    ;;
  *)
    echo "Invalid target: $TARGET"
    echo "Choose one of: check, test, pypi, all"
    exit 1
    ;;
esac

if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required but not found in PATH."
  exit 1
fi

if [[ ! -f "$REPO_ROOT/pyproject.toml" ]]; then
  echo "pyproject.toml not found in $REPO_ROOT"
  exit 1
fi

if [[ ! -x "$VENV_PYTHON" ]]; then
  echo "==> Creating release virtual environment"
  python3 -m venv "$VENV_DIR"
fi

echo "==> Checking git status"
cd "$REPO_ROOT"
git status --short

read -r -p "Proceed with release version ${VERSION}? [y/N] " answer
case "$answer" in
  [yY]|[yY][eE][sS])
    ;;
  *)
    echo "Release cancelled."
    exit 1
    ;;
esac

echo "==> Installing release tools into local venv"
"$VENV_PYTHON" -m pip install --upgrade pip build twine >/dev/null

echo "==> Updating project version"
cd "$REPO_ROOT"
VERSION="$VERSION" "$VENV_PYTHON" - "$REPO_ROOT" <<'PY'
import os
import sys
from pathlib import Path

root = Path(sys.argv[1])
version = os.environ['VERSION']
path = root / 'pyproject.toml'
lines = path.read_text().splitlines()
updated = False
for i, line in enumerate(lines):
    if line.strip().startswith('version = '):
        lines[i] = f'version = "{version}"'
        updated = True
        break
if not updated:
    raise SystemExit('Could not find version line in pyproject.toml')
path.write_text('\n'.join(lines) + '\n')
PY

echo "==> Cleaning previous build artifacts"
rm -rf dist build *.egg-info

echo "==> Building distribution"
"$VENV_PYTHON" -m build

echo "==> Checking package metadata"
"$VENV_PYTHON" -m twine check dist/*

if [[ "$TARGET" == "check" ]]; then
  echo "==> Validation-only mode: build succeeded, no upload was performed."
  exit 0
fi

if [[ "$TARGET" == "test" || "$TARGET" == "all" ]]; then
  echo "==> Uploading to TestPyPI"
  "$VENV_PYTHON" -m twine upload --repository testpypi dist/*
fi

if [[ "$TARGET" == "pypi" || "$TARGET" == "all" ]]; then
  echo "==> Uploading to PyPI"
  "$VENV_PYTHON" -m twine upload dist/*
fi

echo "==> Creating git tag v${VERSION}"
git add pyproject.toml
#git commit -m "Release v${VERSION}"
#git push origin HEAD
#git tag "v${VERSION}"
#git push origin "v${VERSION}"

echo "Release script completed."
echo "Next steps:"
echo "  1. Review the version bump in pyproject.toml"
echo "  2. Commit and push your release branch"
echo "  3. Push the tag if desired"
echo "  4. Create the GitHub release"
