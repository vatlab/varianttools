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
 
# LOG_LEVELS = {'DEBUG': logging.DEBUG,
#               'INFO': logging.INFO,
#               'WARN': logging.WARN,
#               'ERROR': logging.ERROR,
#               'CRITICAL': logging.CRITICAL
#               }
# port = 6000


# logger = logging.getLogger()
# context = zmq.Context()
# socket_fd = context.socket(zmq.SUB)
# socket_fd.bind('tcp://'+os.environ["ZEROMQIP"]+':%i' % port)
# socket_fd.setsockopt(zmq.SUBSCRIBE, b'')
# filehandler = logging.handlers.TimedRotatingFileHandler('log file', 'midnight',1)
# level=logging.DEBUG
# logger.setLevel(level)
# filehandler.setLevel(level)
# formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
# filehandler.setFormatter(formatter)
# logger.addHandler(filehandler)


# while True:
#         topic, message = socket_fd.recv_multipart()
#         pos = topic.find('.')
#         level = topic
#         if pos > 0: level = topic[:pos]
#         if message.endswith('\n'): message = message[:-1]
#         log_msg = getattr(logging, level.lower())
#         if pos > 0: message = topic[pos+1:] + " | " + message
#         log_msg(message)