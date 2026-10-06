# Pinned environment for the time-based lockbox (TIMELOCK.md).
# Build from the pre-registration commit:
#   git checkout <prereg commit>   # git log --diff-filter=A --format=%H -- TIMELOCK.md
#   docker build -t bending-timelock .
#   docker run --rm -v $PWD:/work -w /work bending-timelock python timelock.py gate --stratum A
# The base image is pinned by digest (python:3.11.15-slim as of 2026-08-07, Docker Hub).
# The Docker daemon was not running in the pre-registration environment, so this file was not
# built there; the lockfile versions match that environment exactly.
FROM python:3.11.15-slim@sha256:90744cff8f32887f075c47d747a173ff333e9e98801667af93c357fa9f5e28ff
ENV OMP_NUM_THREADS=1 PYTHONHASHSEED=0
COPY requirements.lock /tmp/requirements.lock
RUN pip install --no-cache-dir -r /tmp/requirements.lock
