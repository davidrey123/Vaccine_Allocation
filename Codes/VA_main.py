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
	S,I,R,D,xall,T = f.optimize_scenario_TS(s,N,P,PM,node2dm,p,same,dm,popdm,S0,I0,bg,cap,1.0,1.0,mode)
	sdm = mode+'_TS_'+e+'_'+s
	f.write_data_dm(sdm,dm,P,S,I,R,D,T)