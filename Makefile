FFTSVDDIR = ../fftsvd

include make-$(HOSTCOMPUTER).inc

OBJS = $(FFTSVDDIR)/Vector.o $(FFTSVDDIR)/ComplexVector.o $(FFTSVDDIR)/Matrix.o \
       $(FFTSVDDIR)/Vector3D.o $(FFTSVDDIR)/Panel.o $(FFTSVDDIR)/Integration.o \
       $(FFTSVDDIR)/VertFace.o $(FFTSVDDIR)/QUI.o $(FFTSVDDIR)/Charge.o \
       $(FFTSVDDIR)/GMRES.o $(FFTSVDDIR)/Cube.o $(FFTSVDDIR)/Tree.o \
       $(FFTSVDDIR)/FFT.o $(FFTSVDDIR)/Preconditioner.o $(FFTSVDDIR)/EquivDensity.o \
       $(FFTSVDDIR)/GreensFunction.o $(FFTSVDDIR)/QuadratureRule.o $(FFTSVDDIR)/LJparameters.o \
       $(FFTSVDDIR)/SurfaceOperator.o $(FFTSVDDIR)/FlatPanel.o $(FFTSVDDIR)/GST.o \
       $(FFTSVDDIR)/TOR.o $(FFTSVDDIR)/FlatIntegration.o $(FFTSVDDIR)/calcpc_GST.o \
       $(FFTSVDDIR)/calcpc_TOR.o $(FFTSVDDIR)/Polynomial.o $(FFTSVDDIR)/FFTSVDpbeAPI.o \
       $(FFTSVDDIR)/QualocationOperator.o $(FFTSVDDIR)/SMatrix.o $(FFTSVDDIR)/SVector.o $(FFTSVDDIR)/ComplexSVector.o \
       PBEproblem.o Unconstrained.o LinConstrained.o BoxConstrained.o Optimizer.o Overlap.o

all:   driverMulti FFTSVDpbeSRF FFTSVDpbeSRF_green getECFDiagonalHessian

obj:  $(OBJS)

clean:
	-rm driverMulti getHessian driverUnc driverLin driverBox FFTSVDpbeSRF getECFDiagonalHessian *.o *.il

fresh: clean all

FFTSVDpbeSRF: FFTSVDpbeSRF.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ FFTSVDpbeSRF.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDpbeSRF_green: FFTSVDpbeSRF_green.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ FFTSVDpbeSRF_green.o $(OBJS) $(LIBDIR) $(LIBS)

getECFDiagonalHessian: getECFDiagonalHessian.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ getECFDiagonalHessian.o $(OBJS) $(LIBDIR) $(LIBS)

getHessian: getHessian.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ getHessian.o $(OBJS) $(LIBDIR) $(LIBS)

getLinear: getLinear.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ getLinear.o $(OBJS) $(LIBDIR) $(LIBS)

driverMulti: driverMulti.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ driverMulti.o $(OBJS) $(LIBDIR) $(LIBS)

driverUnc: driverUnc.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ driverUnc.o $(OBJS) $(LIBDIR) $(LIBS)

driverLin: driverLin.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ driverLin.o $(OBJS) $(LIBDIR) $(LIBS)

driverBox: driverBox.o $(OBJS)
	$(CC) $(LDFLAGS) $(INCLUDE) -o $@ driverBox.o $(OBJS) $(LIBDIR) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c $<
