#!/bin/sh

# copy to Documents folder in my iCloud

if [ $HOME = /Users/junkoda ]; then
  if [ -s cross.pdf ]; then
    cp cross.pdf ~/Library/Mobile\ Documents/com\~apple\~CloudDocs/Documents/
  fi
fi
