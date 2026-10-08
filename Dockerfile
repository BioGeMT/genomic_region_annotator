FROM continuumio/miniconda3:24.7.1-0

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

RUN echo 'precedence ::ffff:0:0/96 100' >> /etc/gai.conf

COPY environment.yml pyproject.toml README.md LICENSE ./
COPY src ./src
COPY scripts ./scripts

RUN conda env update --name base --file environment.yml \
    && conda clean --all --yes \
    && python -m pip install --no-deps . \
    && command -v wget \
    && command -v gunzip

RUN useradd --create-home --shell /usr/sbin/nologin appuser
USER appuser
WORKDIR /work

ENTRYPOINT ["genomic-region-annotator"]
