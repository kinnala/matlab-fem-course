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
            obj.p = p;
            obj.t = t;
            % build affine mappings from (p,t)
            [B, c] = obj.affineMappings(p, t);
            % compute determinants and save them to obj.detB 
            obj.detB = obj.determinant(B);
            % compute inverse transposes of affine mapping
            invBt = obj.inverseTranspose(B, obj.detB);
            % fetch local quadrature rules
            [Q, ~] = obj.localQuadrature();
            % evaluate local basis at quadrature points
            [Phi, Phix, Phiy] = obj.localBasis(Q);
            % map quadrature points to global elements and save to obj.x
            obj.x = obj.F(B, c, Q);
            % transform local basis to global basis and save to obj.phi*
            [obj.phi, obj.phix, obj.phiy] = obj.globalBasis(Phi, Phix, Phiy, invBt);
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
            [~, W] = obj.localQuadrature();
            % initialize
            I = zeros(1, 3*3*size(obj.t,2));
            J = I;
            V = I;
            % loop over local stiffness matrix elements
            for j=1:3
                du = {obj.phix{j}, obj.phiy{j}};
                u = obj.phi{j};
                for i=1:3
                    dv = {obj.phix{i}, obj.phiy{i}};
                    v = obj.phi{i};
                    % save data
                    ix = 3*(i-1)+j;
                    range = (size(obj.t,2)*(ix - 1) + 1):(size(obj.t, 2)*ix);
                    I(range) = obj.t(i,:);
                    J(range) = obj.t(j,:);
                    % weak formulation
                    V(range) = (bilin_form(u,du,v,dv,obj.x)*W').*obj.detB';
                end
            end
            A = sparse(I, J, V, size(obj.p,2), size(obj.p,2));
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
            [~, W] = obj.localQuadrature();
            % initialize
            
            I = zeros(1, 3*size(obj.t,2));
            J = I;
            V = I;
            % loop over local load vector elements
            for i=1:3
                dv = {obj.phix{i}, obj.phiy{i}};
                v = obj.phi{i};
                % save data
                range = (size(obj.t,2)*(i - 1) + 1):(size(obj.t, 2)*i);
                I(range) = obj.t(i,:);
                J(range) = ones(1,size(obj.t,2));
                % weak formulation
                V(range) = (lin_form(v,dv,obj.x)*W').*obj.detB';
            end
            b = full(sparse(I, J, V, size(obj.p,2), 1));
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
            Q = [1/6 2/3 1/6;
                 1/6 1/6 2/3];
            W = [1/6 1/6 1/6];
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
            B = {p(1,t(2,:)) - p(1,t(1,:)), p(1,t(3,:)) - p(1,t(1,:));
                 p(2,t(2,:)) - p(2,t(1,:)), p(2,t(3,:)) - p(2,t(1,:))};
            c = {p(1,t(1,:)), p(2,t(1,:))};
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
            x{1} = bsxfun(@plus, B{1,1}'*Q(1,:) + B{1,2}'*Q(2,:), c{1}');
            x{2} = bsxfun(@plus, B{2,1}'*Q(1,:) + B{2,2}'*Q(2,:), c{2}');
        end
        function detB = determinant(B)
            % Compute the determinants of given 2 x 2 matrices.
            %
            % Input:
            %    B - {2 x 2}(1 x Nel), matrices.
            % Output:
            %    detB - 1 x Nel, determinants of the matrices
            detB = B{1,1}.*B{2,2} - B{1,2}.*B{2,1};
        end
        function invBt = inverseTranspose(B, detB)
            % Compute the inverse transposes of given 2 x 2 matrices.
            %
            % Input:
            %    B - {2 x 2}(1 x Nel), matrices.
            %    detB - (1 x Nel), determinants of the matrices.
            % Output:
            %    invBt - {2 x 2}(1 x Nel), inverse transpose matrices.
            invBt = {B{2,2}./detB,-B{2,1}./detB;
                     -B{1,2}./detB, B{1,1}./detB};
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
            Phi = {1-Q(1,:)-Q(2,:), Q(1,:), Q(2,:)};
            Phix = {-1+0*Q(1,:), 1+0*Q(1,:), 0*Q(1,:)};
            Phiy = {-1+0*Q(1,:), 0*Q(1,:), 1+0*Q(1,:)};
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
            phi = cell(1,length(Phi));
            phix = cell(1,length(Phi));
            phiy = cell(1,length(Phi));
            for k=1:length(Phi)
                phi{k} = repmat(Phi{k}, length(invBt{1,1}), 1);
                phix{k} = invBt{1,1}'*Phix{k} + invBt{1,2}'*Phiy{k};
                phiy{k} = invBt{2,1}'*Phix{k} + invBt{2,2}'*Phiy{k};
            end
        end
    end
end