function[res] = test(N);

    addpath('../mtx/regu');
    addpath('../mtx/SJget');

    % Fix random seed for reproducibility
%    rng(0);

    % Matrices
    matrixNameIdx = 1;
    matrixFunIdx  = 2;
    matrixInfo = { ...
        {'rand'         @(n) rand(n,n) }
        {'vander'       @(n) vander(linspace(0, 1, n)) }
        {'baart'        @baart       }
        {'break1'       @break1      }
        {'break9'       @break9      }
        {'deriv2'       @deriv2      }
        {'devil'        @devil       }
        {'exponential'  @exponential }
        {'foxgood'      @foxgood     }
        {'gravity'      @gravity     }
        {'heat'         @heat        }
        {'phillips'     @phillips    }
        {'random'       @random      }
        {'shaw'         @shaw        }
        {'spikes'       @spikes      }
        {'stewart'      @stewart     }
        {'ursell'       @ursell      } 
        {'wing'         @wing        }
        };

    function[p] = permute_1(A);
        [m,nn] = size(A);
        p = 1:nn;
    end
    
    function[p] = permute_2(A);
        [~, p] = sort(vecnorm(A,2));
    end
    
    function[p] = permute_3(A);
        [~, p] = sort(vecnorm(A,2));
        p = flip(p);
    end
    
    function[p] = permute_4(A);
        [Q,R,p] = qr(A, 'vector');
    end
    
    function[p] = permute_5(A);
        [Q,R,p] = qr(A, 'vector');
        p = flip(p);
    end

    permutes = { ...
        @permute_1
        @permute_2
        @permute_3
        @permute_4
        @permute_5
        };

    
    houses = { ...
        @orth_geqr2
        @orth_golub
        @orth_lapack
        @orth_higham
        };

    function[Q, R, P, deff] = mynone(A);
        [m,nn] = size(A);
        R = A;
        Q = eye(m);
        P = eye(nn);
        deff = zeros(1,nn);
    end
    
    function[Q, R, P, deff] = myqr(A, house);
        [m,nn] = size(A);
        [V,R,T] = householder_qr(A, house);
        Q = eye(m) - V*T*V';
        P = eye(nn);
        deff = zeros(1,nn);
    end
    
    function[Q, R, P, deff] = mypoqr(A, house, def_cri, nrm, nrmA);
        [m,nn] = size(A);
        [V,R,T,deff] = householder_poqr(A, house, def_cri, nrm, nrmA);
        Q = eye(m) - V*T*V';
        P = eye(nn);
    end

    function[Q, R, P, deff] = myrrqr(A, house);
        [m,nn] = size(A);
        [V,R,T,P] = householder_cpqr(A, house);
        Q = eye(m) - V*T*V';
        deff = zeros(1,nn);
    end
    
    function[Q, R, P, deff] = mqr(A);
        [m,nn] = size(A);
        [Q,R] = qr(A);
        P = eye(nn);
        deff = zeros(1,nn);
    end
    
    function[Q, R, P, deff] = mrrqr(A);
        [m,nn] = size(A);
        [Q,R,P] = qr(A);
        deff = zeros(1,nn);
    end
    
    facs = { ...
        @(A, house, def_cri, nrm, nrmA) myqr(A, house)
        @(A, house, def_cri, nrm, nrmA) mypoqr(A, house, def_cri, nrm, nrmA)
        @(A, house, def_cri, nrm, nrmA) myrrqr(A, house)
        @(A, house, def_cri, nrm, nrmA) mqr(A)
        @(A, house, def_cri, nrm, nrmA) mrrqr(A)
        };

    nrms = { ...
        @(A) norm(A, 1)
        @(A) norm(A, 'fro')
        @(A) max(vecnorm(A,2))
        @(A) norm(A, 2)
        @(A) norm(A, 'inf')
        };

    function[cond] = def_cri_1(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < eps * nrmA;
    end

    function[cond] = def_cri_2(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < eps * nrm(A(1:end,1:k));
    end

    function[cond] = def_cri_3(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < eps * nrm(A(1:end,k));
    end

    function[cond] = def_cri_4(R, k, nrm, A, nrmA);
        cond = sum(abs(diag(R(1:k,1:k)))) < eps * nrm(A(1:end,k));
    end

    function[cond] = def_cri_5(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < max(size(A)) * eps(nrmA);
    end

    function[cond] = def_cri_6(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < max(size(A)) * eps(nrm(A(1:end,1:k)));
    end

    function[cond] = def_cri_7(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < max(size(A)) * eps(nrm(A(1:end,k)));
    end

    function[cond] = def_cri_8(R, k, nrm, A, nrmA);
        cond = sum(abs(diag(R(1:k,1:k)))) < max(size(A)) * eps(nrm(A(1:end,k)));
    end

    function[cond] = def_cri_9(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < .1 * eps * nrmA;
    end

    function[cond] = def_cri_10(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < .1 * eps * nrm(A(1:end,1:k));
    end

    function[cond] = def_cri_11(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < .1 * eps * nrm(A(1:end,k));
    end

    function[cond] = def_cri_12(R, k, nrm, A, nrmA);
        cond = sum(abs(diag(R(1:k,1:k)))) < .1 * eps * nrm(A(1:end,k));
    end

    function[cond] = def_cri_13(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < .1 * max(size(A)) * eps(nrmA);
    end

    function[cond] = def_cri_14(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < .1 * max(size(A)) * eps(nrm(A(1:end,1:k)));
    end

    function[cond] = def_cri_15(R, k, nrm, A, nrmA);
        cond = abs( R(k,k) ) < .1 * max(size(A)) * eps(nrm(A(1:end,k)));
    end

    function[cond] = def_cri_16(R, k, nrm, A, nrmA);
        cond = sum(abs(diag(R(1:k,1:k)))) < .1 * max(size(A)) * eps(nrm(A(1:end,k)));
    end

    def_cris = { ...
        @def_cri_1
        @def_cri_2
        @def_cri_3
        @def_cri_4
        @def_cri_5
        @def_cri_6
        @def_cri_7
        @def_cri_8
        @def_cri_9
        @def_cri_10
        @def_cri_11
        @def_cri_12
        @def_cri_13
        @def_cri_14
        @def_cri_15
        @def_cri_16
        };

    function[R2, R_mask] = drop_deff_cols(n, R1, deff1, ignor_deff_cols_i);

        R_mask = 1:n;
        R2 = R1;
        if (sum(deff1) > 0)
            R_mask = R_mask(~deff1);
        else
            if ((ignor_deff_cols_i > 1))
                R_mask = R_mask(diag(R1) > .1 * max(size(R1)) * eps(nrm(R1)));
                R2 = R1(:,R_mask);
            end
        end

    end

    post_facs = { ...
        @(R, house) mynone(R)
        @(R, house) myqr(R, house)
        @(R, house) myrrqr(R, house)
        };

    function[res] = orth_cri1(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        res = size(R3);
        res = res(2);
    end

    function[res] = orth_cri2(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        try
            res = rank(R3);
        catch
            res = -1;
        end
    end

    function[res] = orth_cri3(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        try
            res = cond(R3);
        catch
            res = -1;
        end
    end

    function[res] = orth_cri4(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        res = norm(x - xsol) / norm(xsol);
    end

    function[res] = orth_cri5(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        res = norm(A * x - b);
    end

    function[res] = orth_cri6(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        res = norm(A * x - b) / (norm(A) * norm(x) + norm(b));
    end

    function[res] = orth_cri7(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        res = norm(A' * (A * x - b));
    end

    function[res] = orth_cri8(A, b, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank);
        res = norm(A' * (A * x - b)) / (norm(A).^2);
    end

    orth_cris = { ...
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri1(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri2(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri3(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri4(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri5(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri6(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri7(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        @(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank) orth_cri8(A, rhs, xsol, x, Q, R, P, deff, R2, Q3, R3, P3, deff3, reference_rank),
        };

    % Generate As, xsols and RHSs
    NRHS = 1;
    As = {};
    xsols= rand(N, NRHS);
    RHS = {};
%    e = eps * rand(N, NRHS);
    for m_i = 1 : length(matrixInfo)
        As{m_i} = matrixInfo{m_i}{matrixFunIdx}(N);
        RHS{m_i} = As{m_i} * xsols;% + e;
    end

%    writecell(As, strcat('../exp/', int2str(N), '/As'));
%    writematrix(xsols, strcat('../exp/', int2str(N), '/xsols'));
%    writecell(RHS, strcat('../exp/', int2str(N), '/RHS'));

%    As = readcell(strcat('../exp/', int2str(N), '/As'));
%    xsols = readmatrix(strcat('../exp/', int2str(N), '/xsols'));
%    RHS = readcell(strcat('../exp/', int2str(N), '/RHS'));

    % Main loop
    tic
    k = 0;
    tot = 9444600;
    res = zeros(length(matrixInfo), 5, 4, 5, 16, 5, 2, 3, 8, NRHS);
    for m_i = 1 : 18
        Aorig = As{m_i};
        [m,n] = size(Aorig);
        for p_i = 1 : 1
            permute = permutes{p_i};
            p = permute(Aorig);
            A = Aorig(:,p);
            reference_rank = rank(A);
            for house_i = 4:4 % geqr2, golub, lapack, higham
                house = houses{house_i};
                for fac_i = [2,4,5]%1:5 % qr, poqr, rrqr, mqr, mrrqr
%                    [m_i, p_i, house_i, fac_i, k/tot*100]
                    if (house_i ~= 4 & fac_i > 3) % mqr or mrrqr
                        continue;
                    else
                        fac = facs{fac_i};
                        for def_cri_i = 1:16
                            if (def_cri_i ~= 1 & fac_i ~= 2)
                                continue;
                            else
                                def_cri = def_cris{def_cri_i};
                                for norm_i = 4:4 % 1, fro, 2aprox, 2, inf
                                    if (norm_i ~= 4 & fac_i ~= 2)
                                        continue;
                                    else
                                        nrm = nrms{norm_i};
                                        nrmA = nrm(A);
                                        [Q1, R1, P1, deff1] = fac(A, house, def_cri, nrm, nrmA);
                                        for ignor_deff_cols_i = 1:1 % False, True
                                            [R2, R_mask] = drop_deff_cols(n, R1, deff1, ignor_deff_cols_i);
                                            for post_fac_i = 1:1 % None, mqr, mrrqr
                                                post_fac = post_facs{post_fac_i};
                                                [Q3, R3, P3, deff3] = post_fac(R2, house);
                                                for orthogonality_criteria = 1:3
                                                    orth_cri = orth_cris{orthogonality_criteria};
                                                    res(m_i, p_i, house_i, fac_i, def_cri_i, norm_i, ignor_deff_cols_i, post_fac_i, orthogonality_criteria, 1:NRHS) = orth_cri(A, [], [], [], Q1, R1, P1, deff1, R2, Q3, R3, P3, deff3, reference_rank);
                                                    k = k + 1;
                                                end
                                                for x_i = 1:NRHS
                                                    xsol = xsols(:,x_i);
                                                    rhs = RHS{m_i}(:,x_i);
                                                    x = zeros(n,1);
                                                    y = Q3' * Q1' * rhs;
                                                    if (post_fac_i == 1)
                                                        x(R_mask) = R3 \ y;
                                                    else
                                                        x(R_mask) = R3 \ y(R_mask);
                                                    end
                                                    x(R_mask) = P3 * x(R_mask);
                                                    x = P1 * x;
                                                    for orthogonality_criteria = 4:8
                                                        orth_cri = orth_cris{orthogonality_criteria};
                                                        res(m_i, p_i, house_i, fac_i, def_cri_i, norm_i, ignor_deff_cols_i, post_fac_i, orthogonality_criteria, x_i) = orth_cri(A, rhs, xsol, x, Q1, R1, P1, deff1, R2, Q3, R3, P3, deff3, reference_rank);
    %                                                    (k / tot) * 100.
                                                        k = k + 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
%            writematrix(res(m_i,1:end,1:end,1:end,1:end,1:end,1:end,1:end,1:end,1:end), strcat('../exp/', int2str(N), '/', matrixInfo{m_i}{matrixNameIdx}, '_', int2str(p_i)));
        end
    end
    toc

    fprintf("               matrix                   Forward                                 Backward\n");
    fprintf("                          QR            PAQR           QRCP           QR          PAQR        QRCP\n");
    for m_i = 1:length(matrixInfo)
        fprintf("%20s  ", matrixInfo{m_i}{matrixNameIdx});
        fprintf("%e  ", res(m_i,1,4,4,1,4,1,1,4,1));
        fprintf("%e  ", res(m_i,1,4,2,3,4,1,1,4,1));
        fprintf("%e  ", res(m_i,1,4,5,1,4,1,1,4,1));
        fprintf("%e  ", res(m_i,1,4,4,1,4,1,1,8,1));
        fprintf("%e  ", res(m_i,1,4,2,3,4,1,1,8,1));
        fprintf("%e\n", res(m_i,1,4,5,1,4,1,1,8,1));
    end 
end
