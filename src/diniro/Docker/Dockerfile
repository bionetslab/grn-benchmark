FROM python:3.9
RUN apt-get update
RUN pip install --upgrade pip
ENV PYTHONUNBUFFERED 1
RUN mkdir /code2
WORKDIR /code2
COPY requirements.txt /code2/
RUN mkdir /umap-master
COPY umap-master/ /umap-master
WORKDIR /umap-master
RUN  pip install -r requirements.txt
RUN python setup.py install
WORKDIR /code2
RUN pip install -r requirements.txt
RUN pip install django-crispy-forms
CMD python manage.py makemigrations
CMD python manage.py migrate
CMD python manage.py collectstatic
COPY . /code2/
