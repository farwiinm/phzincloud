# ── Stage 1: Builder ──────────────────────────────────────────────────────────
FROM python:3.11-slim-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        libffi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir --target=/app/deps -r requirements.txt

# ── Stage 2: Runtime ──────────────────────────────────────────────────────────
FROM python:3.11-slim-bookworm AS runtime

LABEL maintainer="Fathima Farwin Mohamed Milhan"
LABEL project="pH-ZinCloud"

WORKDIR /app

COPY --from=builder /app/deps /app/deps
ENV PYTHONPATH="/app/deps"

COPY src/ ./src/
COPY tests/ ./tests/
COPY validate_parser.py .
COPY batch_parse.py .

RUN mkdir -p results logs data/raw/test_proteins data/reference

ENV PDB_INPUT_DIR="data/raw/test_proteins"
ENV OUTPUT_CSV="results/results.csv"
ENV CUTOFF="5.0"
ENV LOG_LEVEL="INFO"
ENV PYTHONUNBUFFERED="1"

CMD ["python", "src/pipeline.py"]
