FROM python:3.9

LABEL version="1.0.2" maintainer="Nicolas Vannieuwkerke <nicolas.vannieuwkerke@ugent.be>"

RUN apt update && apt install tabix python3-pip -y

RUN pip install setuptools==57.5.0

COPY ./requirements.txt .
RUN pip install -r ./requirements.txt

CMD ["python"]
