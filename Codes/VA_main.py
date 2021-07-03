import VA_functions as f

N,P,PM,node2dm = f.read_data_node('data_node.txt')
p,same = f.read_data_mobility('data_mobility.txt')
dm,popdm,S0,I0,bg,cap = f.read_data_dm('data_dm.txt')

#---without budget sharing
mode = 'none'

#---with budget sharing
#mode = 'flow2cap'

for ind in range(100):
	s = 'X104_20_'+str(ind+1)+'.txt'
	
	#---solve VA using Thompson Sampling algorithm
	S,I,R,D,xall,T = f.optimize_scenario_TS(s,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,1.0,1.0,mode)
	sdm = mode+'_TS_'+s
	f.write_data_dm(sdm,dm,P,S,I,R,D,T)
	
	#---solve VA using Population-Based algorithm
	S,I,R,D,xall,T = f.optimize_scenario_PB(s,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,1.0,1.0,mode)
	sdm = mode+'_PB_'+s
	f.write_data_dm(sdm,dm,P,S,I,R,D,T)	
	
	#---solve VA using Moving-Average algorithm
	S,I,R,D,xall,T = f.optimize_scenario_MA(s,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,1.0,1.0,mode)
	sdm = mode+'_MA_'+s
	f.write_data_dm(sdm,dm,P,S,I,R,D,T)		
	
	#---solve VA using GreedY algorithm
	S,I,R,D,xall,T = f.optimize_scenario_GY(s,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,1.0,1.0,mode)
	sdm = mode+'_GY_'+s
	f.write_data_dm(sdm,dm,P,S,I,R,D,T)		
