TSCFLAGS=-t esnext -m esnext

all: cloth.js

cloth.js: _cloth.js
	gsed -e '/^import/d' _cloth.js > cloth.js

%.js: %.ts
	-tsc $(TSCFLAGS) $<

.PHONY: clean all

clean:
	-rm _*.js
	-rm cloth.js
