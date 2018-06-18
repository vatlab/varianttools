from celery import Celery


app = Celery(
    'task_receiver',
    broker='amqp://guest@localhost//',
    backend='rpc://',
    include=['variant_tools.celery_main.task_receiver'],
    accept_content="pickle",
    task_serializer="pickle",
    result_serializer="pickle"
)

app.conf.update(
    task_serializer='pickle',
    accept_content=['pickle'],  # Ignore other content
    result_serializer='pickle'
)
