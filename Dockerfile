FROM python:3.14-slim@sha256:9aaa474097d094040649f2162a627c0049b408f52138d980916788d80f04802f

WORKDIR /workspace
COPY requirements.txt .
RUN python -m pip install --no-cache-dir -r requirements.txt && python -m pip check
COPY notebooks/fibril_reproducibility.ipynb notebooks/

CMD ["jupyter", "execute", "notebooks/fibril_reproducibility.ipynb", "--timeout=60"]
