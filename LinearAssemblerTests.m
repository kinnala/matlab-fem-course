classdef LinearAssemblerTests < matlab.unittest.TestCase
    % Write "runtests LinearAssemblerTests" to run.
    methods (Test)
        function testStaticLocalQuadrature(test)
            p = @(x,y) 7*x.^2-5*x.*y-6*y.^2+2*x-10*y+55;
            [Q,W] = LinearAssembler.localQuadrature();
            test.verifyEqual(size(Q,1),2);
            test.verifyEqual(sum(p(Q(1,:),Q(2,:)).*W), 625/24, ...
                             'AbsTol', 1e-5)
        end
        function testStaticAffineMappings(test)
            % single triangle
            p = [0 0; 1 0; 0 1]';
            t = [1 2 3]';
            [B, c] = LinearAssembler.affineMappings(p, t);
            test.verifyEqual(size(B),[2 2]);
            test.verifyEqual(size(c),[1 2]);
            test.verifyEqual([B{1,1} B{1,2}; B{2,1} B{2,2}],eye(2));
            test.verifyEqual([c{1,1} c{1,2}],[0 0]);
            % two triangles
            midpts={[1/3 1/3], [2/3 2/3]};
            p = [0 0; 1 0; 0 1; 1 1]';
            t = [1 2 3; 2 3 4]';
            [B, c] = LinearAssembler.affineMappings(p, t);
            for itr=1:2
                x=[1/3 1/3]';
                Bk=[B{1,1}(itr) B{1,2}(itr); B{2,1}(itr) B{2,2}(itr)];
                ck=[c{1,1}(itr) c{1,2}(itr)]';
                test.verifyEqual(inv(Bk)*Bk,eye(2));
                test.verifyEqual(Bk*x+ck,midpts{itr}','AbsTol', 1e-8);
            end
        end
        function testStaticF(test)
            % check that output size is correct
            p = [0 0; 1 0; 0 1]';
            B = {[1 1], [0 0]; [0 0], [1 1]};
            c = {[0 1], [1 0]};
            x = LinearAssembler.F(B, c, p);
            test.verifyEqual(size(x{1}),[2 3]);
            test.verifyEqual(size(x{2}),[2 3]);
            % identity mapping
            p = [0 0; 1 0; 0 1]';
            B = {1, 0; 0, 1};
            c = {0, 0};
            x = LinearAssembler.F(B, c, p);
            test.verifyEqual(p(1,:),x{1});
            test.verifyEqual(p(2,:),x{2});
            % scaling
            p = [0 0; 1 0; 0 1]';
            B = {2, 0; 0, 2};
            c = {0, 0};
            x = LinearAssembler.F(B, c, p);
            test.verifyEqual(2*p(1,:),x{1});
            test.verifyEqual(2*p(2,:),x{2});
            % translate
            p = [0 0; 1 0; 0 1]';
            B = {1, 0; 0, 1};
            c = {1, 1};
            x = LinearAssembler.F(B, c, p);
            test.verifyEqual(p(1,:)+1,x{1});
            test.verifyEqual(p(2,:)+1,x{2});
        end
        function testStaticDeterminant(test)
            B = {[1 1], [0 0]; [0 0], [1 1]};
            detB = LinearAssembler.determinant(B);
            test.verifyEqual(length(detB),2);
            test.verifyEqual(detB,[1 1]);
        end
        function testStaticInverseTranspose(test)
            B = {10, 2.6; 1.2, 25.1};
            invBt = LinearAssembler.inverseTranspose(B,...
                LinearAssembler.determinant(B));
            test.verifyEqual(inv(cell2mat(B))',cell2mat(invBt), ...
                'AbsTol', 1e-8);
        end
        function testStaticLocalBasis(test)
            p = [0 0; 1 0; 0 1]';
            [phi,phix,phiy]=LinearAssembler.localBasis(p);
            % check sizes
            test.verifyEqual(length(phi),3);
            test.verifyEqual(length(phix),3);
            test.verifyEqual(length(phiy),3);
            % test that values at nodes sum to one
            test.verifyEqual(sum(phi{1}),1);
            test.verifyEqual(sum(phi{2}),1);
            test.verifyEqual(sum(phi{3}),1);
            % check derivatives approximately
            eps = 1e-10;
            point = 0.3333333;
            p = [point+eps point; point point+eps; point point]';
            [phi,phix,phiy]=LinearAssembler.localBasis(p);
            for j=1:3
                test.verifyEqual((phi{j}(1)-phi{j}(3))/eps,phix{j}(1),...
                    'AbsTol',1e-4);
                test.verifyEqual((phi{j}(2)-phi{j}(3))/eps,phiy{j}(1),...
                    'AbsTol',1e-4);
            end
        end
        function testStaticGlobalBasis(test)
            % test out sizes
            B = {[1 1], [0 0]; [0 0], [1 1]};
            invBt = LinearAssembler.inverseTranspose(B,...
                LinearAssembler.determinant(B));
            p = [0.333 0.333; 0 1]';
            [Phi,Phix,Phiy]=LinearAssembler.localBasis(p);
            [phi, phix, phiy] = ...
                LinearAssembler.globalBasis(Phi, Phix, Phiy, invBt);
            test.verifyEqual(length(phi),3);
            for j=1:3
                test.verifyEqual(size(phi{j}),[2 2]);
                test.verifyEqual(size(phix{j}),[2 2]);
                test.verifyEqual(size(phiy{j}),[2 2]);
            end
            % one element, one point
            B = {1, 0; 0, 1};
            invBt = LinearAssembler.inverseTranspose(B,...
                LinearAssembler.determinant(B));
            p = [0.333 0.333]';
            [Phi,Phix,Phiy]=LinearAssembler.localBasis(p);
            [phi, phix, phiy] = ...
                LinearAssembler.globalBasis(Phi, Phix, Phiy, invBt);
            test.verifyEqual(phi,Phi);
            test.verifyEqual(phix,Phix);
            test.verifyEqual(phiy,Phiy);
        end
        function testConstructor(test)
            % just test that everything is filled
            l = LinearAssembler([0,0;1,0;0,1;1,1]',[1,2,3;2,3,4]');
            test.verifyFalse(isempty(l.x));
            test.verifyFalse(isempty(l.p));
            test.verifyFalse(isempty(l.t));
            test.verifyFalse(isempty(l.detB));
            test.verifyFalse(isempty(l.phi));
            test.verifyFalse(isempty(l.phix));
            test.verifyFalse(isempty(l.phiy));
        end
        function testAssemblerBilinear(test)
            % check that the area of (0,0) (1,0) (0,1) integrates to 0.5
            l = LinearAssembler([0,0;1,0;0,1]',[1,2,3]');
            M = l.assembleBilinear(@(u,du,v,dv,x) u.*v);
            test.verifyEqual(sum(sum(full(M))),0.5,'AbsTol',1e-10);
            
            % check that the function x in square (0,1)^2 integrates to 0.5
            p = [0,0;1,0;1,1;0,1]';
            t = [1,2,4;2,3,4]';
            l = LinearAssembler(p,t);
            M = l.assembleBilinear(@(u,du,v,dv,x) x{1}.*u.*v);
            test.verifyEqual(sum(sum(full(M))),0.5,'AbsTol',1e-10);
            
            % check that the function y in square (0,1)^2 integrates to 0.5
            M = l.assembleBilinear(@(u,du,v,dv,x) x{2}.*u.*v);
            test.verifyEqual(sum(sum(full(M))),0.5,'AbsTol',1e-10);
            
            % check that the function x^2+x*y+y^2+x+y in square (0,1)^2
            % integrates to 23/12
            M = l.assembleBilinear(@(u,du,v,dv,x) (x{1}.^2+x{1}.*x{2}+x{2}.^2+x{1}+x{2}).*u.*v);
            test.verifyEqual(sum(sum(full(M))),23/12,'AbsTol',1e-10);
            
            % check derivatives
            
            % derivative of constant is zero
            M = l.assembleBilinear(@(u,du,v,dv,x) du{1}.*v);
            test.verifyEqual(ones(1,4)*M*ones(4,1),0,'AbsTol',1e-10);
            M = l.assembleBilinear(@(u,du,v,dv,x) du{2}.*v);
            test.verifyEqual(ones(1,4)*M*ones(4,1),0,'AbsTol',1e-10);
            
            % derivative of x (or y) with respect to x (or y) should be 1
            M = l.assembleBilinear(@(u,du,v,dv,x) du{1}.*v);
            test.verifyEqual(ones(1,4)*M*p(1,:)',1,'AbsTol',1e-10);
            M = l.assembleBilinear(@(u,du,v,dv,x) du{2}.*v);
            test.verifyEqual(ones(1,4)*M*p(2,:)',1,'AbsTol',1e-10);
        end
        function testAssemblerLinear(test)
            p = [0,0;1,0;1,1;0,1]';
            t = [1,2,4;2,3,4]';
            l = LinearAssembler(p,t);
            f = l.assembleLinear(@(v,dv,x) x{1}.*v);
            test.verifyEqual(size(f),[4 1]);
            test.verifyEqual(sum(f),0.5,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) x{2}.*v);
            test.verifyEqual(sum(f),0.5,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) dv{1});
            test.verifyEqual(f'*ones(4,1),0,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) dv{2});
            test.verifyEqual(f'*ones(4,1),0,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) dv{2});
            test.verifyEqual(f'*p(1,:)',0,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) dv{2});
            test.verifyEqual(f'*p(2,:)',1,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) dv{1});
            test.verifyEqual(f'*p(2,:)',0,'AbsTol',1e-10);
            f = l.assembleLinear(@(v,dv,x) dv{1});
            test.verifyEqual(f'*p(1,:)',1,'AbsTol',1e-10);
        end
    end
end