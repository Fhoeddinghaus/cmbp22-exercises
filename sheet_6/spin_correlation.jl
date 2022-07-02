function spin_correlation(N,ψ, n,m; H_sc = false)
    if n == m
        # algorithm for hamiltonian not working for (n,n) coupling, manually:
        return 3/4 * ψ' * ψ # 1/4 for each α=x,y,z
    end
    if H_sc == false
        # calculate the Hamilton matrix for the coupling
        Js_sc = reset_Js(N)
        # only coupling between (n,m) for all α
        Js_sc[1][n,m] = -1
        Js_sc[2][n,m] = -1
        Js_sc[3][n,m] = -1
        H_sc = calculate_hamilton_matrix(Js_sc, N)
    end # else: H_sc is given, e.g.for multiple calculations for the same (n,m) for different states
    return (ψ' * H_sc * ψ)[1] # 1×1 matrix to scalar
end