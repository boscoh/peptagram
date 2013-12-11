

text = open('style.css').read()

import scss
print scss.Scss(scss_opts={'compress':True}).compile(scss_file='style.css')

# from scss import parser
# print parser.parse(text)