from celery import Celery
import os

server=os.environ["RABBITMQIP"]

app = Celery(
    'task_receiver',
    # broker='amqp://guest@localhost//',
    broker='amqp://guest@'+os.environ["RABBITMQIP"]+"//",
    backend='rpc://',
    # backend='redis://guest@'+os.environ["RABBITMQIP"]+':6379/0',
    include=['variant_tools.celery_main.task_receiver'],
    accept_content="pickle",
    task_serializer="pickle",
    result_serializer="pickle"
    
)
app.conf.broker_heartbeat = 0

app.conf.update(
    task_serializer='pickle',
    accept_content=['pickle'],  # Ignore other content
    result_serializer='pickle'
)
