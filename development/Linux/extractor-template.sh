#!/usr/bin/env bash
#
# Varianttools Installer
# Gao Wang and Bo Peng Copyright 2013
# 
ARCH=x64
BUNDLE=
EXE=vtools

# VTOOLS_KEEP_TEMP=1

if [ -z "$BASH" ]; then
   # $- expands to the current options so things like -x get passed through
   if [ ! -z "$-" ]; then
      opts="-$-"
   fi

   # dash flips out of $opts is quoted, so don't.
   exec /usr/bin/env bash $opts "$0" "$@"
   echo "Unable to restart with bash shell"
   exit 1
fi

set -e

# Determines whether the user land is 32 or 64-bit.
get_arch() {
   # First byte is the ELF magic number.  The 5th byte is whether it's
   # a 32 or 64-bit machine (1 or 2, respectively).  See `man elf` for
   # details.
   local ELF_MAGIC=7f

   if [ "`od -N1 -An -t x1 < /bin/sh | tr -d ' '`" != "$ELF_MAGIC" ]; then
      exit 1
   fi

   local arch=`od -j4 -N1 -An -t u1 < /bin/sh | tr -d ' '`

   case $arch in
      1)
         echo "x86"
	 exit 0
	 ;;
      2)
         echo "x64"
         exit 0
	 ;;
      *)
         exit 1
         ;;
   esac
}

# Determines if path is relative.
# 0 if relative, otherwise 1.
is_relative() {
    local path="$1"
    shift

    [ "${path:0:1}" != "/" ]
    return
}

# Determines if a program is in the user's PATH.
# 0 if found, else 1
internal_which() {
   local binary="$1"

   for dir in `echo $PATH | tr ":" "\n"`; do
      if [ -s "$dir/$binary" -a -x "$dir/$binary" ]; then
         return 0
      fi
   done

   return 1
}

# create tmpdir VTOOLS_TEMP and extract bundle there
extract() {
   VTOOLS_TEMP=`mktemp -d /tmp/vtools.XXXXXX`

   if [ ! -d "$VTOOLS_TEMP" ]; then
      echo "Unable to create temporary directory."
      exit 1
   fi
   
   echo -e "Extracting files ...\n"

   (cd $VTOOLS_TEMP && sed -e '1,/^exit$/d' "$1" | tar xzf -)
}

bootstrap() {
   bash $1/INSTALL $1
   ret=$?
   if [ $ret != 0 ]; then
      exit $ret
   fi
   return 0
}

# Called when the script exits
# VTOOLS_TEMP removed unless VTOOLS_KEEP_TEMP is set
on_exit() {
   if [ -e "$VTOOLS_TEMP" -a -z "$VTOOLS_KEEP_TEMP" ]; then
      rm -rf "$VTOOLS_TEMP"
   fi
}

# Execute when exit signal is received
trap on_exit EXIT

main() {
   local fullpath=$0

   if [ "`get_arch`" != $ARCH ]; then
      echo "This program is for Linux $ARCH and does not match "
      echo "that of the current architecture."
      exit 1
   fi

   if is_relative $fullpath; then
      fullpath=$PWD/$fullpath
   fi
   extract $fullpath
   bootstrap $VTOOLS_TEMP/$BUNDLE
   if internal_which $EXE; then
       echo "Installation complete!"
   else
       echo "[WARNING] $EXE is not installed to your path! "
       echo "You will not be able to execute command $EXE. "
       exit 1
   fi
}

main $@
exit
