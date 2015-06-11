INC = -I./include
LIBDIR = lib
SRCDIR = src
BUILDDIR = build
LIB = -lxdrfile -larmadillo -lboost_program_options
DEBUG = -g
LFLAGS += $(DEBUG) -Wall
CFLAGS += $(DEBUG) -Wall -c
LIBFLAGS = -c -static -fpic
TARGET = bin/z_sfg

SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

$(TARGET): $(OBJECTS)
	@echo $(SOURCES)
	@echo " Linking..."
	@echo " $(CXX) $^ -o $(TARGET) $(LIB)"; $(CXX) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

