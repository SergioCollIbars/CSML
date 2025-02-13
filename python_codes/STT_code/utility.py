# This script is holds all miscellaneous functions that reshape or restructure python objects

def reshape_STT(X,order,stateSize): 
    X_state = X[0:stateSize] # Extract state
    phi_out_list = []
    
    if order >= 1: # Extract STM
        second_start = stateSize + stateSize**2
        phi_first = X[stateSize:second_start].reshape((stateSize,stateSize))
        phi_out_list.append(phi_first)
        
    if order >= 2: # Extract 2nd Order STT
        third_start = second_start + stateSize**3
        phi_second = X[second_start:third_start].reshape((stateSize,stateSize,stateSize))
        phi_out_list.append(phi_second)
        
    if order >= 3:
        fourth_start = third_start + stateSize**4
        phi_third = X[third_start:fourth_start].reshape((stateSize,stateSize,stateSize,stateSize))
        phi_out_list.append(phi_third)
        
    if order >= 4:
        phi_fourth = X[fourth_start:(fourth_start+stateSize**5)].reshape((stateSize,stateSize,stateSize,stateSize,stateSize))
        phi_out_list.append(phi_fourth)
        
    return X_state,phi_out_list