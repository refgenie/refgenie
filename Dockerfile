# Multi-stage: the image must be self-contained and reproducible from a bare
# checkout. `deploy-api.yml` runs `docker build` from a plain checkout where
# refgenie/server/webui/ does not exist (it is gitignored, built by CI); left
# single-stage, that would silently deploy a UI-less API.
FROM node:26-slim AS web
WORKDIR /build
# package.json/package-lock.json first, then the rest of the source, so the
# `npm ci` layer stays cached across frontend source edits.
COPY frontend/package.json frontend/package-lock.json ./frontend/
RUN cd frontend && npm ci
COPY frontend ./frontend
COPY pyproject.toml ./
# Writes /build/refgenie/server/webui/ (index.html, hashed _app/*, build-info.json).
RUN cd frontend && npm run build

FROM python:3.14-slim
LABEL authors="Nathan Sheffield, Michal Stolarczyk, Oleksandr Khoroshevskyi"

RUN apt-get update && apt-get install -y \
    gcc \
    libpq-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app
# The build stage is the only source of truth for the web UI bundle -- see
# .dockerignore, which excludes refgenie/server/webui from the `COPY . /app`
# context above so a maintainer's local build can never leak into this image.
COPY --from=web /build/refgenie/server/webui /app/refgenie/server/webui

RUN pip install --upgrade pip && pip install uv
RUN uv pip install ".[server]" --system

COPY deployment/config.yaml /config.yaml
ENV REFGENIE_DB_CONFIG_PATH=/config.yaml

EXPOSE 80
CMD ["uvicorn", "refgenie.server.main:app", "--host", "0.0.0.0", "--port", "80"]
