# Releasing scTap

This document gives the exact release workflow for publishing the package to PyPI.

<i>Use this exact flow for a new public release.</i>

### 1. Confirm metadata

Check that [pyproject.toml](pyproject.toml) contains:

```toml
[project]
name = "scTap"
version = "0.0.2"
```

Also confirm the README install instructions match:

```bash
pip install scTap
```

### 2. Validate the build locally

From the project root:

```bash
./release.sh 0.0.2 check
```

This validates the package build without uploading anything.

### 3. Commit the release state

```bash
git add .
git commit -m "Prepare scTap 0.0.2"
git push origin main
```

### 4. Create the tag that triggers the GitHub workflow

```bash
git tag v0.0.2
git push origin v0.0.2
```

The workflow in [.github/workflows/publish.yml](.github/workflows/publish.yml) is configured to run on tags matching `v*`.

### 5. Wait for GitHub Actions to publish

The workflow will:

- set up Python
- build the package
- publish to PyPI

### 6. Verify the release

After the workflow completes, confirm that:

- the package appears on PyPI as `scTap`
- version `0.0.2` is visible
- installation works with:

```bash
pip install scTap
```

