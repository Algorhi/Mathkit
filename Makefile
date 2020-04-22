#-------------------------------------------------------------------------------
# Makefile for the mathkit library
#-------------------------------------------------------------------------------


include ../reentry.mk

INCLUDES = -lm

C = $(CC) $(CFLAGS) $(INCLUDES) -Iinclude

Source := $(wildcard src/*.c)
OBJS := $(patsubst %.c,%.o,$(Source))
.PHONY:all
all: libmk.a 

#-------------------------------------------------------------------------------
# the mathkit library:
#-------------------------------------------------------------------------------
.PHONY:demo
demo:test
test:test.c libmk.a
	$(CC) $(CFLAGS) -o $@ $^ -Iinclude $(LDFLAGS)
	@echo mathkit demo successfully built. Type ./test to run demo problem.


libmk.a: $(OBJS)
	$(ARCHIVE) $@ $^
	-$(RANLIB) $@
	@echo libmk.a successfully built.
src/%.o:src/%.c 
	$(C) -c $< -o $@
.PHONY:purge clean
purge:
	- $(RM) libmk.a $(OBJS) test
clean:
	- $(RM) $(CLEAN) libmk.a
