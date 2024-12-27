# Compiler
FC = mpif90

# Compiler flags
FFLAGS = -ffree-line-length-none -O3

# Source files
SRCS = utils.f90 main.f90

# Executable name
EXEC = genic

# Default target
all: $(EXEC)

# Rule to build the executable
$(EXEC): $(SRCS)
	$(FC) $(FFLAGS) $(SRCS) -o $(EXEC)

# Clean target to remove the executable
clean:
	rm -f $(EXEC) && rm *.mod
