#!/usr/bin/env python

# asciify.py -- converts bitmap image to ASCII art -- Sep 9, 2010 -- Georgy Samsonidze
# based on ASCII Art Maker by Steven Kay
# http://stevendkay.wordpress.com/2009/09/08/generating-ascii-art-from-photographs-in-python/
# the ASCII grayscale table is taken from IMG2ASCII by Ueli Weiss
# http://img2ascii.sourceforge.net/

def main(argv = None):
  if argv is None:
    argv = sys.argv
  argc = len(argv)
  self = "asciify.py"
  if argv[0][-len(self):] != self:
    sys.exit("\n   Rename script to %s\n" % self)
  if argc != 8:
    sys.exit("\n   Usage: %s mode fimage fproc awidth aoffset charw charh\n" % self +
             "      mode = 1 - preprocess & write | 2 - read & display | 3 - normal\n" +
             "    fimage = name of file to read image from (format supported by PIL)\n" +
             "     fproc = name of file to write/read preprocessed image to/from\n" +
             "    awidth = width of asciified image in characters\n" +
             "   aoffset = number of blank characters in front of asciified image\n" +
             "     charw = width of a character in pixels (for aspect ratio)\n" +
             "     charh = height of a character in pixels (for aspect ratio)\n")

  try:
    mode = int(argv[1])
  except:
    mode = -1
  fimage = argv[2]
  fproc = argv[3]
  try:
    awidth = int(argv[4])
  except:
    awidth = -1
  try:
    aoffset = int(argv[5])
  except:
    aoffset = -1
  try:
    charw = int(argv[6])
  except:
    charw = -1
  try:
    charh = int(argv[7])
  except:
    charh = -1

  if mode < 1 or mode > 3:
    sys.exit("\n   Error: illegal value %s for mode\n" % argv[1])
  if awidth < 1:
    sys.exit("\n   Error: illegal value %s for awidth\n" % argv[4])
  if aoffset < 0:
    sys.exit("\n   Error: illegal value %s for aoffset\n" % argv[5])
  if charw < 1:
    sys.exit("\n   Error: illegal value %s for charw\n" % argv[6])
  if charh < 1:
    sys.exit("\n   Error: illegal value %s for charh\n" % argv[7])

  if mode == 1 or mode == 3:
    try:
      from PIL import Image
    except:
      sys.exit("\n   Error: unable to load module PIL\n")
    try:
      img = Image.open(fimage)
      (oldwidth, oldheight) = img.size
      newwidth = awidth
      newheight = int(float(oldheight * newwidth * charw) / float(oldwidth * charh) + 0.5)
      img = img.resize((newwidth, newheight), Image.ANTIALIAS)
      img = img.convert('L')
      nx = img.size[0]
      ny = img.size[1]
      dat = []
      for iy in range(ny):
        dat.append([])
        for ix in range(nx):
          dat[iy].append(img.getpixel((ix, iy)))
    except:
      sys.exit("\n   Error: unable to read file %s\n" % fimage)

  if mode == 1:
    try:
      h = open(fproc, 'w')
      s = ' %i %i\n' % (nx, ny)
      h.write(s)
      for iy in range(ny):
        s = ''
        for ix in range(nx):
          s = s + ' %i'  % dat[iy][ix]
        s = s + '\n'
        h.write(s)
      h.close()
    except:
      sys.exit("\n   Error: unable to write file %s\n" % fproc)

  if mode == 2:
    try:
      h = open(fproc, 'r')
      s = h.readline()
      t = s.split()
      nx = int(t[0])
      ny = int(t[1])
      dat = []
      for iy in range(ny):
        s = h.readline()
        t = s.split()
        dat.append([])
        for ix in range(nx):
          dat[iy].append(int(t[ix]))
      h.close()
    except:
      sys.exit("\n   Error: unable to read file %s\n" % fproc)

  if mode == 2 or mode == 3:
    import random
    greyscale = [ ' ',
                  ' ',
                  '`',
                  '`',
                  '\'.',
                  '\'.',
                  '!,-^~',
                  '!,-^~',
                  ':_',
                  ';{}',
                  '()/<>\\',
                  '*?',
                  '"=',
                  '%+il',
                  '17,t',
                  '3Ccjr',
                  '&2IJovxz',
                  '$Lu',
                  '6OTYZ',
                  '45Genswy',
                  '0ADFPVf',
                  '9SUahp',
                  '8Xbdk',
                  'QRq',
                  'BEKWgm',
                  '@',
                  '@',
                  'H',
                  '#',
                  'N',
                  'M' ]
    greylength = len(greyscale)
    s = ''
    for iy in range(ny):
      for j in range(aoffset):
        s = s + ' '
      for ix in range(nx):
        lum = int(float((255 - dat[iy][ix]) * (greylength - 1)) / float(255) + 0.5)
        possibles = greyscale[lum]
        s = s + possibles[random.randint(0, len(possibles) - 1)]
      s = s + '\n'
    print(s)

  return 0

if __name__ == "__main__":
   import sys
   sys.exit(main())

