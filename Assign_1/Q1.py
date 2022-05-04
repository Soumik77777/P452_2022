import library


matA = [[1,-1,4,0,2,9],[0,5,-2,7,8,4],[1,0,5,7,3,-2],[6,-1,2,3,0,8],[-4,2,0,5,-5,3],[0,7,-1,5,4,2]]
vec = [19,2,13,-7,-9,2]

print("Solution using Gauss Jordan Method:")
soln1 = library.gauss_jordan(matA, vec)
print("[a1, a2, a3, a4, a5, a6]= ", soln1)
print("")
print("Solution using LU decomposition Method:")
soln2 = library.LUSolve(matA, vec)
print("[a1, a2, a3, a4, a5, a6]= ", soln2)



'''
Solution using Gauss Jordan Method:
[a1, a2, a3, a4, a5, a6]=  [-1.8042964518987756, 0.7454404218772226, 3.98178356910756, -1.5696545989650006, 2.0707049157011412, 0.16457697399371995]

Solution using LU decomposition Method:
[a1, a2, a3, a4, a5, a6]=  [-1.8042964518987756, 0.7454404218772226, 3.98178356910756, -1.5696545989650006, 2.0707049157011412, 0.16457697399371995]

'''

