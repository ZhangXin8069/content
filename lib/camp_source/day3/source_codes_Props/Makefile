TOPDIR=/dssg/home/acct-phyww/phyww/chroma/chroma_single_gcc-8.5/
CHROMA=$(TOPDIR)/chroma_single
CONFIG=$(CHROMA)/bin/chroma-config


MGLIBS= -lwilsonmg -lqopqdp -lqdp_common -lqdp_int -lqdp_f -lqdp_d -lqdp_df -lqdp_f2 -lqdp_d2 -lqdp_df2 -lqdp_f3 -lqdp_d3 -lqdp_df3 -lqdp_fn -lqdp_dn -lqdp_dfn -lqio -llime -lqla_c99 -lqla_cmath -lqla_d2 -lqla_d3 -lqla_d -lqla_df2 -lqla_df3 -lqla_df -lqla_dfn -lqla_dn -lqla_dq2 -lqla_dq3 -lqla_dq -lqla_dqn -lqla_f2 -lqla_f3 -lqla_f -lqla_fn -lqla_int -lqla_q2 -lqla_q3 -lqla_q -lqla_qn -lqla_random
MGLDFLAGS=-L$(TOPDIR)/qla/lib -L$(TOPDIR)/lib -L$(TOPDIR)/qopqdp/lib
MGCXXFLAGS=-L$(TOPDIR)/qla/lib -L$(TOPDIR)/qdp/lib -L$(TOPDIR)/qopqdp/lib

CXX=$(shell $(CONFIG) --cxx) -D__HIP_PLATFORM_HCC__
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I$(TOPDIR)/qla/include -I. $(MGCXXFLAGS) 
#CXXFLAGS=-D_PDF_ $(shell $(CONFIG) --cxxflags) -I$(TOPDIR)/qla/include -I. $(MGCXXFLAGS) 
LDFLAGS=$(shell $(CONFIG) --ldflags) $(MGLDFLAGS)
LIBS=$(shell $(CONFIG) --libs) $(MGLIBS)

HDRS=grid_source.h \
     inline_qpropadd_w.h \
     inline_seqsource_fast_w.h \
     simple_baryon_seqsrc_w.h \

OBJS=chroma.o \
     grid_source.o \
     inline_qpropadd_w.o \
     inline_seqsource_fast_w.o \
     simple_baryon_seqsrc_w.o \

chroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

hmc: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)


%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -rf chroma hmc $(OBJS) *~
