import ROOT

matrix = ROOT.TMatrixD(2,2)
matrix[0][0]=1
matrix[0][1]=2
matrix[1][0]=2
matrix[1][1]=3

matrix.Print()




eigen = ROOT.TMatrixDEigen(matrix)

eigen_value_matrix = eigen.GetEigenValues()
eigen_value_matrix.Print()

eigen_vector_matrix = eigen.GetEigenVectors()
eigen_vector_matrix.Print()

for i in range(2):
    print 'Lambda '+str(i+1)+': '+str(eigen_value_matrix[i][i])


