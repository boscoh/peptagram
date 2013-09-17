from __future__ import print_function
import logging

log_level = 30
logging.basicConfig(level=log_level)
print('base', log_level)

# root_logger = logging.getLogger('')
# root_logger.setLevel(5)
# print('ROOT', root_logger.level)

logger = logging.getLogger('example')
print('example', logger.root.level)
logger.info('haha')