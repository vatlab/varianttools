import zmq
import logging
import logging.handlers
import os

port = 6001
ctx = zmq.Context()
sub = ctx.socket(zmq.SUB)
sub.bind('tcp://*:%i' % port)
sub.setsockopt(zmq.SUBSCRIBE, b"")
level=logging.DEBUG
logging.basicConfig(level=level)


while True:
    level, message = sub.recv_multipart()
    message = message.decode('ascii')
    if message.endswith('\n'):
        # trim trailing newline, which will get appended again
        message = message[:-1]
    log = getattr(logging, level.lower().decode('ascii'))
    log(message)
