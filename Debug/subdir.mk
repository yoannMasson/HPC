################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../message.c \
../mpi_hello_world.c 

CPP_SRCS += \
../Analytic.cpp \
../CrankNicolson.cpp \
../FTCS.cpp \
../Laasonen.cpp \
../Solver.cpp \
../Test.cpp \
../matrix.cpp \
../vector.cpp 

OBJS += \
./Analytic.o \
./CrankNicolson.o \
./FTCS.o \
./Laasonen.o \
./Solver.o \
./Test.o \
./matrix.o \
./message.o \
./mpi_hello_world.o \
./vector.o 

C_DEPS += \
./message.d \
./mpi_hello_world.d 

CPP_DEPS += \
./Analytic.d \
./CrankNicolson.d \
./FTCS.d \
./Laasonen.d \
./Solver.d \
./Test.d \
./matrix.d \
./vector.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


