#!/usr/bin/env bash
set -euo pipefail

julia --project=. benchmark/generate_report.jl
