#!/bin/sh
cd "$(dirname "$0")/.."
exec python3 tests/check.py
