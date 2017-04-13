classdef LinearAssembler
    % LinearAssembler Finite element assembly using P1 triangles
    %
    %    Provided a mesh (p,t), this class performs FEM assembly, that is,
    %    the class transforms bilinear forms into sparse matrices.
    %
    %    An example for assembling the mass matrix M:
    %
    %        a = LinearAssembler(p,t);
    %        M = a.assembleBilinear(@(u,du,v,dv,x) u.*v);
    %
    %    In the documentation we use the following notation:
    %
    %       Nve                - number of VErtices in the mesh
    %       Nel                - number of ELements in the mesh
    %       Nqp                - number of Quadrature Points
    %
    %       2 x Nve            - matrix with 2 rows and Nve columns
    %       Nel x Nqp          - matrix with Nel rows and Nqp columns
    %       {1 x 3}(Nel x Nqp) - cell matrix with 1 row and 3 columns
    %                            with elements of type Nel x Nqp
    %
    %    See <a href="matlab:doc('LinearAssembler')">doc LinearAssembler</a> for more information.
    properties
        p    % 2 x Nve, mesh vertices 
        t    % 3 x Nel, mesh elements
        detB % 1 x Nel, determinants of affine mappings in each element
        phi  % {1 x 3}(Nel x Nqp), basis at global quadrature points
        phix % {1 x 3}(Nel x Nqp), basis x derivatives at -||-
        phiy % {1 x 3}(Nel x Nqp), basis y derivatives at -||-
        x    % {1 x 2}(Nel x Nqp), locations of global quadrature points
    end
    methods
        function obj = LinearAssembler(p, t)
            % Fill in all properties of the class given the mesh (p,t).
            %
            % Input:
            %    p - 2 x Nve, mesh vertices
            %    t - 3 x Nel, mesh elements
            % Output:
            %    obj - the constructed assembler object
            %
        end
        function A = assembleBilinear(obj, bilin_form)
            % Transform a bilinear form @(u,du,v,dv,x) into a matrix.
            %
            % Input:
            %    bilin_form - a function handle describing the bilinear
            %                 form and that has the following parameters:
            %                 @(u,du,v,dv,x).
            % Output:
            %    A - Nve x Nve, a sparse matrix corresponding to the
            %                   bilinear form.
            %
        end
        function b = assembleLinear(obj, lin_form)
            % Transform a linear form @(v,dv,x) into a matrix.
            %
            % Input:
            %    lin_form - a function handle describing the linear
            %               form and that has the following parameters:
            %               @(v,dv,x).
            % Output:
            %    b - Nve x 1, a vector corresponding to the linear form.
            %
        end
    end
    methods (Static)
        function [Q, W] = localQuadrature()
            % Return second-order accurate quadrature.
            %
            % Input:
            %    None
            % Output:
            %    Q - 2 x Nqp, each column is a quadrature point on the 
            %                 reference element.
            %    W - 1 x Nqp, each value is a quadrature weight of the
            %                 corresponding quadrature point.
        end
        function [B, c] = affineMappings(p, t)
            % Build affine mappings F(X)=BX+c in all elements.
            %
            % Input:
            %    p - 2 x Nve, mesh vertices
            %    t - 3 x Nel, mesh elements
            % Output:
            %    B - {2 x 2}(1 x Nel), affine mapping matrices in every
            %                          element of the mesh (p,t).
            %    c - {1 x 2}(1 x Nel), affine mapping constants in every
            %                          element.
        end
        function x = F(B, c, Q)
            % Map points defined on the reference element affinely.
            %
            % Input:
            %    B - {2 x 2}(1 x Nel), affine mapping matrices.
            %    c - {1 x 2}(1 x Nel), affine mapping constants.
            %    Q - 2 x Nqp, each column is a point defined on the 
            %                 reference element.
            % Output:
            %    x - {1 x 2}(Nel x Nqp), points Q mapped through the mapping
            %                            F(Q)=BQ+c.
        end
        function detB = determinant(B)
            % Compute the determinants of given 2 x 2 matrices.
            %
            % Input:
            %    B - {2 x 2}(1 x Nel), matrices.
            % Output:
            %    detB - 1 x Nel, determinants of the matrices
        end
        function invBt = inverseTranspose(B, detB)
            % Compute the inverse transposes of given 2 x 2 matrices.
            %
            % Input:
            %    B - {2 x 2}(1 x Nel), matrices.
            %    detB - (1 x Nel), determinants of the matrices.
            % Output:
            %    invBt - {2 x 2}(1 x Nel), inverse transpose matrices.
        end
        function [Phi, Phix, Phiy] = localBasis(Q)
            % Evaluate local basis functions and their derivatives.
            %
            % Input:
            %    Q - 2 x Nqp, local points inside the reference triangle.
            %
            % Output:
            %    Phi - {1 x 3}(1 x Nqp), basis function values
            %    Phix - {1 x 3)(1 x Nqp), basis function x derivatives
            %    Phiy - {1 x 3)(1 x Nqp), basis function y derivatives
        end
        function [phi, phix, phiy] = globalBasis(Phi, Phix, Phiy, invBt)
            % Transform local basis into global basis.
            %
            % Input:
            %    Phi - {1 x 3}(1 x Nqp), local basis function values
            %    Phix - {1 x 3)(1 x Nqp), local basis function x derivatives
            %    Phiy - {1 x 3)(1 x Nqp), local basis function y derivatives
            %    invBt - {2 x 2}(1 x Nqp), inverse transpose matrices.
            %
            % Output:
            %    phi - {1 x 3}(Nel x Nqp), global basis function values
            %    phix - {1 x 3}(Nel x Nqp), global basis function
            %                               x derivatives
            %    phiy - {1 x 3}(Nel x Nqp), global basis function 
            %                               y derivatives
        end
    end
end
