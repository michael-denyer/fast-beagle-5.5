CC ?= cc
CFLAGS ?= -O2 -g
# Part of correctness, not tuning: see the plan's Architecture section.
override CFLAGS += -std=c11 -Wall -Wextra -Werror -ffp-contract=off -fno-fast-math -fwrapv -Isrc -Ithird_party/libdeflate
# Honours the `#pragma omp simd` lines on the HMM value loops. It links no OpenMP runtime.
override CFLAGS += -fopenmp-simd
# fdlibm is kept as upstream wrote it; these warnings flag its style, not bugs.
FDLIBM_CFLAGS := -Wno-dangling-else -Wno-sign-compare
# libdeflate 1.25 as plink2 vendors it (bgen=plink2 must compress as plink2 does), built with its own flags.
LIBDEFLATE_CFLAGS := -O2
JAVA ?= java
JAVAC ?= javac
PREFIX ?= /usr/local
LDLIBS += -lm -pthread
# htslib from Homebrew on macOS; on Linux the system package needs no flags.
# HTSLIB_PREFIX=<dir> uses the htslib in <dir> instead, as the conda build does.
HTSLIB_PREFIX ?= $(shell brew --prefix htslib 2>/dev/null)
ifneq ($(HTSLIB_PREFIX),)
override CFLAGS += -I$(HTSLIB_PREFIX)/include
LDFLAGS += -L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
endif
LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lhts $(LDLIBS)

JCOMPAT_OBJ := $(patsubst src/%.c,build/obj/%.o,$(wildcard src/jcompat/*.c src/jcompat/fdlibm/*.c))
JCOMPAT_FIXTURES := random math numbers utf8 parse parseint pqueue search
LIBDEFLATE_OBJ := $(patsubst %.c,build/obj/%.o,$(wildcard third_party/libdeflate/lib/*.c third_party/libdeflate/lib/*/*.c))
BEAGLE_OBJ := $(sort $(patsubst src/%.c,build/obj/%.o,$(wildcard src/*/*.c)) $(JCOMPAT_OBJ) $(LIBDEFLATE_OBJ))

.PHONY: all install check-jcompat check-bgen-unit check-records check-bgen-files check-vcf-index check-tbi check-tracker check-oom check-interval check-markers check-block-reader check-sliding-window check-snv-perms check-piece-size java-trace clean
.SECONDARY:
.DELETE_ON_ERROR:
all: build/beagle

build/beagle: $(BEAGLE_OBJ)
	$(LINK)

# Installs build/beagle as fast-beagle, since Beagle's own packages install beagle.
install: build/beagle
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 build/beagle $(DESTDIR)$(PREFIX)/bin/fast-beagle

check-jcompat: $(JCOMPAT_FIXTURES:%=build/jcompat/%.diff)

build/obj/jcompat/fdlibm/%.o: src/jcompat/fdlibm/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(FDLIBM_CFLAGS) -MMD -MP -c -o $@ $<

build/obj/third_party/%.o: third_party/%.c
	@mkdir -p $(@D)
	$(CC) $(LIBDEFLATE_CFLAGS) -MMD -MP -c -o $@ $<

build/obj/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -MMD -MP -c -o $@ $<

-include $(BEAGLE_OBJ:.o=.d) build/obj-piece1/imp/imputed_writer.d

# build/beagle with one reference marker per imputation work item, in its own
# object directory so that build/beagle keeps the default.
build/obj-piece1/imp/imputed_writer.o: src/imp/imputed_writer.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -DPIECE_RECORDS=1 -MMD -MP -c -o $@ $<

build/beagle-piece1: $(filter-out build/obj/imp/imputed_writer.o,$(BEAGLE_OBJ)) build/obj-piece1/imp/imputed_writer.o
	$(LINK)

check-piece-size: build/beagle build/beagle-piece1
	tests/check-piece-size.sh build/beagle build/beagle-piece1

build/jcompat/%_fixture: tests/jcompat/%_fixture.c $(JCOMPAT_OBJ)
	@mkdir -p $(@D)
	$(LINK)

build/jcompat/pqueue_fixture: build/obj/beagleutil/comp_hap_queue.o build/obj/blbutil/utilities.o

check-bgen-unit: build/bgen/quantise_test build/bgen/info_test build/bgen/pack_test
	./build/bgen/quantise_test
	./build/bgen/info_test
	./build/bgen/pack_test

build/bgen/%_test: tests/bgen/%_test.c $(filter-out build/obj/main/main.o,$(BEAGLE_OBJ))
	@mkdir -p $(@D)
	$(LINK)

check-records: build/output/record_fixture
	python3 -B tests/check_records.py

build/output/record_fixture: tests/output/record_fixture.c $(filter-out build/obj/main/main.o,$(BEAGLE_OBJ))
	@mkdir -p $(@D)
	$(LINK)

check-bgen-files: build/output/bgen_files_fixture
	python3 -B tests/check_bgen_files.py

build/output/bgen_files_fixture: tests/output/bgen_files_fixture.c build/obj/bgen/bgen_files.o build/obj/main/run_outputs.o build/obj/blbutil/utilities.o $(JCOMPAT_OBJ)
	@mkdir -p $(@D)
	$(LINK)

check-vcf-index: build/output/vcf_index_test
	@mkdir -p build/output/vcf-index
	./build/output/vcf_index_test build/output/vcf-index

build/output/vcf_index_test: tests/output/vcf_index_test.c build/obj/main/vcf_index.o build/obj/blbutil/str_set.o \
        build/obj/blbutil/utilities.o $(JCOMPAT_OBJ)
	@mkdir -p $(@D)
	$(LINK)

check-tbi: build/beagle build/output/vcf_index_test
	tests/check-tbi.sh

check-tracker: build/beagleutil/tracker_test
	./build/beagleutil/tracker_test

build/beagleutil/tracker_test: tests/beagleutil/tracker_test.c build/obj/beagleutil/comp_hap_queue.o \
        build/obj/blbutil/utilities.o $(JCOMPAT_OBJ)
	@mkdir -p $(@D)
	$(LINK)

check-oom: build/blbutil/oom_test
	./build/blbutil/oom_test

build/blbutil/oom_test: tests/blbutil/oom_test.c build/obj/blbutil/utilities.o $(JCOMPAT_OBJ)
	@mkdir -p $(@D)
	$(LINK)

check-interval: build/vcf/interval_it_test
	./build/vcf/interval_it_test

check-markers: build/vcf/markers_test
	./build/vcf/markers_test

build/vcf/markers_test: tests/vcf/markers_test.c $(filter-out build/obj/main/main.o,$(BEAGLE_OBJ))
	@mkdir -p $(@D)
	$(LINK)

check-block-reader: build/vcf/block_reader_test
	./build/vcf/block_reader_test

check-sliding-window: build/vcf/sliding_window_test
	./build/vcf/sliding_window_test

build/vcf/sliding_window_test: tests/vcf/sliding_window_test.c $(filter-out build/obj/main/main.o,$(BEAGLE_OBJ))
	@mkdir -p $(@D)
	$(LINK)

check-snv-perms: build/vcf/snv_perms_test
	./build/vcf/snv_perms_test

# The scheduling test includes the implementation to interpose unlock.
# It links the engine without main.o, since the test defines its own main.
build/vcf/block_reader_test: tests/vcf/block_reader_test.c src/vcf/block_reader.c $(filter-out build/obj/main/main.o build/obj/vcf/block_reader.o,$(BEAGLE_OBJ))
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(filter-out src/vcf/block_reader.c,$^) -lhts $(LDLIBS)

# The table test includes the implementation to read its static table.
build/vcf/snv_perms_test: tests/vcf/snv_perms_test.c src/vcf/marker.c $(filter-out build/obj/main/main.o build/obj/vcf/marker.o,$(BEAGLE_OBJ))
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(filter-out src/vcf/marker.c,$^) -lhts $(LDLIBS)

build/vcf/interval_it_test: tests/vcf/interval_it_test.c build/obj/vcf/interval_it.o \
        build/obj/beagleutil/chrom_interval.o build/obj/beagleutil/chrom_ids.o build/obj/blbutil/str_set.o \
        build/obj/blbutil/utilities.o $(JCOMPAT_OBJ)
	@mkdir -p $(@D)
	$(LINK)

build/jcompat/classes/JcompatFixtures.class: tests/jcompat/JcompatFixtures.java
	@mkdir -p $(@D)
	$(JAVAC) -d $(@D) $<

build/jcompat/%.java.txt: build/jcompat/classes/JcompatFixtures.class
	$(JAVA) -cp $(<D) JcompatFixtures $* > $@

# Each C fixture must print the Java output exactly. Fixtures with inputs read them from it.
# The diff prints to the log; the .diff file is a stamp written only on a match.
build/jcompat/%.diff: build/jcompat/%_fixture build/jcompat/%.java.txt
	@rm -f $@
	./$< < build/jcompat/$*.java.txt > build/jcompat/$*.c.txt
	diff -u build/jcompat/$*.java.txt build/jcompat/$*.c.txt
	@touch $@
	@echo "PASS jcompat $*"

# The Java source with java/trace.patch applied: dumps trace seams under -Dbeagle.trace=<dir>.
java-trace: build/java-trace/classes/main/Main.class

build/java-trace/classes/main/Main.class: java/trace.patch $(shell find java/src -name '*.java')
	rm -rf build/java-trace
	mkdir -p build/java-trace
	cp -R java/src build/java-trace/src
	patch -d build/java-trace/src -p1 < java/trace.patch
	$(JAVAC) -nowarn -d build/java-trace/classes $$(find build/java-trace/src -name '*.java')

clean:
	rm -rf build
