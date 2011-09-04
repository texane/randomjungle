#!/usr/bin/env sh

RJ_DIR=$HOME/repo/randomjungle/trunk ;
RJ_O="" ;

for f in \
$RJ_DIR/src/lr/*.c \
; do
cpp_name=$f;
o_name=${cpp_name/.c/.o};
gcc \
-Wall \
-Wno-unused-variable \
-Wno-unused-but-set-variable \
-Wno-unused-function \
-Wno-sign-compare \
-Wno-char-subscripts \
-Wunknown-pragmas \
-Wwrite-strings \
-I/usr/include/libxml2 \
-I$RJ_DIR/src/library \
-I$RJ_DIR/src/lr \
-I$RJ_DIR/src \
-c -o $o_name $cpp_name ;
RJ_O="$RJ_O $o_name" ;
done

for f in \
$RJ_DIR/src/library/*.cpp \
; do
cpp_name=$f;
o_name=${cpp_name/.cpp/.o};
g++ \
-Wall \
-Wno-unused-variable \
-Wno-unused-but-set-variable \
-Wno-unused-function \
-Wno-sign-compare \
-Wno-char-subscripts \
-Wunknown-pragmas \
-I/usr/include/libxml2 \
-I$RJ_DIR/src/library \
-I$RJ_DIR/src/lr \
-I$RJ_DIR/src \
-c -o $o_name $cpp_name ;
RJ_O="$RJ_O $o_name" ;
done

rm librj.a ;
ar q librj.a $RJ_O ;
ranlib librj.a ;
