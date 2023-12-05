KUBECONFIG=config-rtis
CONTEXT= <your context here>
KUBECTL=/usr/bin/kubectl --context=${CONTEXT}
build: 
	echo "Container is built and published by github actions. https://github.com/tmolteno/tart_cal"

delete: 
	-${KUBECTL} delete -f signal_cal_cronjob.yaml  -f rhodes_cal_cronjob.yaml
	${KUBECTL} delete -f data-volume-persistentvolumeclaim.yaml

create: 
	${KUBECTL} create -f data-volume-persistentvolumeclaim.yaml
	${KUBECTL} get persistentvolumeclaim
	
cron:
	-${KUBECTL} delete -f signal_cal_cronjob.yaml -f rhodes_cal_cronjob.yaml
	${KUBECTL} create -f rhodes_cal_cronjob.yaml
	${KUBECTL} create -f signal_cal_cronjob.yaml 
	${KUBECTL} get cronjob
	
status:
	${KUBECTL} get pods
	${KUBECTL} get jobs
	${KUBECTL} get persistentvolumeclaim
	${KUBECTL} get cronjob

trigger:
	${KUBECTL} create job --from=cronjob/signal-cal signal-triggered

.PHONY: k9s

k9s:
	 ./k9s_u/k9s --context=${CONTEXT}
