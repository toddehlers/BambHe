#
# Makefile for creating bambHe executable
#

CC=ifort
CLFAGS=-O3
MAIN_OBJS=bambHe.o ball_main.o hexcylinder_main.o hexcylinder2_main.o hexcylinder3_main.o hexdisc_main.o rounddisc_main.o
MOD_OBJS=ball_module.o hexcylinder_module.o hexcylinder2_module.o hexcylinder3_module.o hexdisc_module.o rounddisc_module.o

bambHe_linux: $(MAIN_OBJS) $(MOD_OBJS)
	$(CC) -o bambHe_linux $(MAIN_OBJS) $(MOD_OBJS)
	cp bambHe_linux ../

bambHe.o: bambHe.f90
	$(CC) $(CFLAGS) -c bambHe.f90

ball_main.o: ball_main.f ball_module.o
	$(CC) $(CFLAGS) -c ball_main.f
ball_module.o: ball_module.f
	$(CC) $(CFLAGS) -c ball_module.f

hexcylinder_main.o: hexcylinder_main.f hexcylinder_module.o
	$(CC) $(CFLAGS) -c hexcylinder_main.f
hexcylinder_module.o: hexcylinder_module.f
	$(CC) $(CFLAGS) -c hexcylinder_module.f

hexcylinder2_main.o: hexcylinder2_main.f hexcylinder2_module.o
	$(CC) $(CFLAGS) -c hexcylinder2_main.f
hexcylinder2_module.o: hexcylinder2_module.f
	$(CC) $(CFLAGS) -c hexcylinder2_module.f

hexcylinder3_main.o: hexcylinder3_main.f hexcylinder3_module.o
	$(CC) $(CFLAGS) -c hexcylinder3_main.f
hexcylinder3_module.o: hexcylinder3_module.f
	$(CC) $(CFLAGS) -c hexcylinder3_module.f

hexdisc_main.o: hexdisc_main.f hexdisc_module.o
	$(CC) $(CFLAGS) -c hexdisc_main.f
hexdisc_module.o: hexdisc_module.f
	$(CC) $(CFLAGS) -c hexdisc_module.f

rounddisc_main.o: rounddisc_main.f rounddisc_module.o
	$(CC) $(CFLAGS) -c rounddisc_main.f
rounddisc_module.o: rounddisc_module.f
	$(CC) $(CFLAGS) -c rounddisc_module.f


clean:
	rm -f *.mod
	rm -f *.o