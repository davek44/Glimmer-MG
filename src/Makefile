TGT= $@


all:
	@ TGT=objs
	@ $(dosubdirs)
	@ TGT=libs
	@ $(dosubdirs)
	@ TGT=progs
	@ $(dosubdirs)


install: all
	@ $(dosubdirs)

clean:
	@ $(dosubdirs)


tester: 
	@ $(dosubdirs)


regression: 
	@ $(dosubdirs)

LOCAL_WORK = $(shell cd ..; pwd)
# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.glm
