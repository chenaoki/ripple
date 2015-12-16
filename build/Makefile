CC=/usr/bin/mpicxx

CFLAGS=-std=c++11 -O3 \
	-I../src \
	-isystem /home/tomii/simenv/lifev \
	-isystem /home/tomii/simenv/lifev/build \
	-isystem /home/tomii/simenv/trilinos/install/include \
	-isystem /usr/include/suitesparse \
	-isystem /usr/include/boost \
	`pkg-config --cflags opencv`

LFLAGS=-std=c++11 -O3

LIFEV_LIB_DIR=/home/tomii/simenv/lifev/build/lifev/

LIBS=-L/home/tomii/simenv/trilinos/install/lib -rdynamic $(LIFEV_LIB_DIR)/heart/liblifevheart.a $(LIFEV_LIB_DIR)/core/liblifevcore.a -L/home/tomii/simenv/trilinos/install/lib -lml -lifpack -lamesos -lanasazitpetra -lModeLaplace -lanasaziepetra -lanasazi -lthyraepetra -lthyracore -lrtop -lbelostpetra -lbelosepetra -lbelos -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetrakernels -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -lkokkosalgorithms -lkokkoscontainers -laztecoo -lzoltan -lepetraext -ltriutils -lepetra -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lkokkoscore -lumfpack -lamd -lhdf5 -lz -lparmetis -lmetis -llapack -lblas -lparmetis -lmetis -llapack -lblas -Wl,-rpath,/home/tomii/simenv/trilinos/install/lib `pkg-config --libs opencv`

OBJS:=main.o \
	Heart.o \
	HeartDiffusionFunctor.o \
	HeartDiffusionModel.o \
	HeartIonicData.o

ripple: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o ripple $(LIBS)

%.o : ../src/%.cpp
	$(CC) $(CFLAGS) $(VERSION) -c $< -o $@

clean:
	rm *.o ripple
