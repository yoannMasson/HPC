################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../labs/PIComputation.c \
../labs/circular.c \
../labs/dotProduct.c \
../labs/hello.c \
../labs/message.c \
../labs/messageFromTask.c \
../labs/messageSendAndRecv.c \
../labs/pingPong.c 

OBJS += \
./labs/PIComputation.o \
./labs/circular.o \
./labs/dotProduct.o \
./labs/hello.o \
./labs/message.o \
./labs/messageFromTask.o \
./labs/messageSendAndRecv.o \
./labs/pingPong.o 

C_DEPS += \
./labs/PIComputation.d \
./labs/circular.d \
./labs/dotProduct.d \
./labs/hello.d \
./labs/message.d \
./labs/messageFromTask.d \
./labs/messageSendAndRecv.d \
./labs/pingPong.d 


# Each subdirectory must supply rules for building sources it contributes
labs/%.o: ../labs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


