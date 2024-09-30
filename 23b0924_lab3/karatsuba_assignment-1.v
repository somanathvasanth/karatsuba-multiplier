/* 32-bit simple karatsuba multiplier */

/*32-bit Karatsuba multipliction using a single 16-bit module*/

module iterative_karatsuba_32_16(clk, rst, enable, A, B, C);
    input clk;
    input rst;
    input [31:0] A;
    input [31:0] B;
    output [63:0] C;
    
    input enable;
    
    
    wire [1:0] sel_x;
    wire [1:0] sel_y;
    
    wire [1:0] sel_z;
    wire [1:0] sel_T;
    
    
    wire done;
    wire en_z;
    wire en_T;
    
    
    wire [32:0] h1;
    wire [32:0] h2;
    wire [63:0] g1;
    wire [63:0] g2;
    
    assign C = g2;
    reg_with_enable #(.N(64)) Z(.clk(clk), .rst(rst), .en(en_z), .X(g1), .O(g2) );  // Fill in the proper size of the register
    reg_with_enable #(.N(33)) T(.clk(clk), .rst(rst), .en(en_T), .X(h1), .O(h2) );  // Fill in the proper size of the register
    
    iterative_karatsuba_datapath dp(.clk(clk), .rst(rst), .X(A), .Y(B), .Z(g2), .T(h2), .sel_x(sel_x), .sel_y(sel_y), .sel_z(sel_z), .sel_T(sel_T), .en_z(en_z), .en_T(en_T), .done(done), .W1(g1), .W2(h1));
    iterative_karatsuba_control control(.clk(clk),.rst(rst), .enable(enable), .sel_x(sel_x), .sel_y(sel_y), .sel_z(sel_z), .sel_T(sel_T), .en_z(en_z), .en_T(en_T), .done(done));
    
endmodule

module iterative_karatsuba_datapath(clk, rst, X, Y, T, Z, sel_x, sel_y, en_z, sel_z, en_T, sel_T, done, W1, W2);
    input clk;
    input rst;
    input [31:0] X;    // input X
    input [31:0] Y;    // Input Y
    input [32:0] T;    // input which sums X_h*Y_h and X_l*Y_l (its also a feedback through the register)
    input [63:0] Z;    // input which calculates the final outcome (its also a feedback through the register)
    output  [63:0] W1;  // Signals going to the registers as input
    output  [32:0] W2;  // signals hoing to the registers as input
    

    input [1:0] sel_x;  // control signal 
    input [1:0] sel_y;  // control signal 
    
    input en_z;         // control signal 
    input [1:0] sel_z;  // control signal 
    input en_T;         // control signal 
    input [1:0] sel_T;  // control signal 
    
    input done;   
          // Final done signal
          wire[15:0] X_l;
          wire [15:0] Y_l;
          wire [15:0] X_h;
          wire [15:0] Y_h;
    assign  X_l=X[15:0];
    assign  Y_l=Y[15:0];
    assign  X_h=X[31:16];
    assign  Y_h=Y[31:16];

    wire[16:0] X_sum;
    wire [16:0] Y_sum;
    adder_Nbit #(.N(16)) xsum (.a(X_h), .b(X_l), .cin(1'b0), .S(X_sum[15:0]), .cout(X_sum[16]));
    adder_Nbit #(.N(16)) ysum (.a(Y_h), .b(Y_l), .cin(1'b0), .S(Y_sum[15:0]), .cout(Y_sum[16]));
    reg [15:0] x_to_multiply;
reg [15:0] y_to_multiply;
reg [31:0] xshifted;
reg [31:0] yshifted;

    
    
    reg [32:0]Wreg;
    reg [63:0]W1reg;
    reg [63:0]W1reg1;
    reg [32:0]W2reg;
    reg [63:0]W2reg1;
    wire [31:0]W;
    wire [33:0] Wnow;
    reg [33:0] Wnowreg;
    reg [32:0] Wnowreg1;
       
        reg [33:0] onlyshifted;
    wire [32:0] sumint1;
 reg [32:0] yshifted1;
wire [33:0] sumint2;
wire ccc;
reg [33:0] yshifted11;
reg [33:0] newx;
 wire ddd;
reg [33:0] xshifted11;
reg [33:0] newy; 
 wire eee;
reg [63:0] shift_w_by_32bit;
reg [63:0] shift_t_by_16bits;
 wire x;
wire xx;
wire pp;
reg [32:0] neww;
wire xxx;
       wire yyy;
reg [33:0]newT;
wire [33:0] answer;
mult_16 multiply (
    .X(x_to_multiply),
    .Y(y_to_multiply),
    .Z(W[31:0])
);



adder_Nbit #(.N(32)) shiftsum (.a(W[31:0]), .b(xshifted[31:0]), .cin(1'b0), .S(sumint1[31:0]), .cout(sumint1[32]));
adder_Nbit #(.N(33)) shiftsum11 (.a(sumint1[32:0]), .b(yshifted1[32:0]), .cin(1'b0), .S(sumint2[32:0]), .cout(sumint2[33]));
 
 //wnow 3 unnay 1 cheyali
 reg [33:0] wnow1st;
 reg [33:0] wnow2nd;
 
adder_Nbit #(.N(34)) shiftsum1111111 (.a(wnow1st), .b(wnow2nd), .cin(1'b0), .S(Wnow[33:0]), .cout(eee)); 
reg [63:0]w11st;
reg [63:0]w12nd;

 adder_Nbit #(.N(64)) wsum2 (.a(w11st), .b(w12nd), .cin(1'b0), .S(W1[63:0]), .cout(xx));
 reg [32:0]w21st;
 reg [32:0]w22nd;
adder_Nbit #(.N(33)) wsum212 (.a(w21st), .b(w22nd), .cin(1'b0), .S(W2[32:0]), .cout(pp));

subtract_Nbit #(.N(34)) asdfff (.a(Wnow[33:0]), .b({1'b0,T[32:0]}), .cin(1'b0), .S(answer[33:0]), .ov(xxx), .cout_sub(yyy)); 






     always@(*) begin
    
    
    case (sel_x)
    2'b01: begin
    x_to_multiply=X_l[15:0];
    end
    2'b10 :begin
     x_to_multiply=X_h[15:0];
    end
    2'b11: begin
     x_to_multiply=X_sum[15:0];
    end
    default:begin
    end
    endcase
    case (sel_y)
    2'b01: begin
     y_to_multiply=Y_l[15:0];
    end
    2'b10: begin
    y_to_multiply=Y_h[15:0];
    end
    2'b11: begin
     y_to_multiply=Y_sum[15:0];
    end
     default:begin
    end
    endcase




    case (X_sum[16])
    1'b1:begin
        case(Y_sum[16])
        1'b1:
        begin
         onlyshifted = {2'b01, 32'b0};

         xshifted = {X_sum[15:0], 16'b0};

         yshifted = {Y_sum[15:0], 16'b0};

        
        
       
         yshifted1={1'b0,xshifted[31:0]};
        
        wnow1st=sumint2[33:0];
        wnow2nd=onlyshifted[33:0];
        
       
        end
        1'b0:
        begin
            
             yshifted11={2'b0,Y_sum[15:0],16'b0};
            
             newx={18'b0,X_sum[15:0]};
           wnow1st=newx[33:0];
           wnow2nd=yshifted11[33:0];
           
          
        end
        endcase
        end
    1'b0:begin
          case(Y_sum[16])
        1'b1:
        begin
             
             xshifted11={2'b00,X_sum[15:0],16'b0000000000000000};
            
             newy = {18'b000000000000000000, Y_sum[15:0]};
           wnow1st=newy[33:0];
           wnow2nd=xshifted11[33:0];
           
        end
        1'b0:
        begin
             Wnowreg1[33:0]={2'b0,W[31:0]};
             wnow1st=34'b0;
             wnow2nd=Wnowreg1;
        end
        endcase
    end
    
    endcase
    



    case (en_z)
    1'b1:begin
        case(sel_z)
    2'b01:begin
     W1reg1[31:0]=W[31:0];
     W1reg1[63:32]={32{1'b0}};
     w11st=W1reg1;
     w12nd=64'b0;
    end
    2'b10:begin
    
     // Concatenates W[31:0] with 32 bits of zeros
    w11st=Z[63:0];
    w12nd={W[31:0], 32'b0};
   
    
    end
    2'b11:begin
    
     shift_t_by_16bits={15'b0,T[32:0],16'b0};
    w11st=Z[63:0];
    w12nd={15'b0,T[32:0],16'b0};
   
    end
        default:begin
            w11st=63'b0;
            w12nd=63'b0;

    end
        endcase
    end
    default: begin
    end
    endcase



    case (en_T)
        1'b1:begin
    case (sel_T)
    2'b01:begin
     W2reg1[32:0]={1'b0,W[31:0]};
     w21st=33'b0;
     w22nd=W2reg1;
    end
     2'b10:begin
        
         neww={1'b0,W[31:0]};
    w21st=T[32:0];
    w22nd=neww[32:0];
    end
     2'b11:begin
       
       
        newT={1'b0,T[32:0]};
       
    
     W2reg1[32:0]=answer[32:0];
     w21st=33'b0;
     w22nd=W2reg1;    
    end
         default:begin
            w21st=33'b0;
            w22nd=33'b0;

    end
        endcase
    end
    default:begin
        
    end
    endcase



     end
    
   
    
    //-------------------------------------------------------------------------------------------------
    
    // Write your datapath here
    //--------------------------------------------------------

endmodule


module iterative_karatsuba_control(clk,rst, enable, sel_x, sel_y, sel_z, sel_T, en_z, en_T, done);
    input clk;
    input rst;
    input enable;
    
    output reg [1:0] sel_x;
    output reg [1:0] sel_y;
    
    output reg [1:0] sel_z;
    output reg [1:0] sel_T;    
    
    output reg en_z;
    output reg en_T;
    
    
    output reg done;
    
    reg [5:0] state, nxt_state;
    parameter S0 = 6'b000001;
    parameter S1 = 6'b000010;
    parameter S2 = 6'b000100;
    parameter S3 = 6'b001000;
    parameter S4 = 6'b010000;
    parameter S5 = 6'b100000;

    
    
      // initial state
   // <define the rest of the states here>

    always @(posedge clk) begin
        if (rst) begin
            state <= S0;
        end
        else if (enable) begin
            state <= nxt_state;
        end
    end
    

    always@(*) begin
        case(state) 
            S0: 
                begin
				sel_x=2'b00;
                sel_y=2'b00;
                sel_z=2'b00;
                sel_T=2'b00;
                en_T=1'b1;
                en_z=1'b1;
                done=1'b0;
                nxt_state=S1;
                	// Write your output and next state equations here
                end
            S1 :
                  begin
				sel_x=2'b01;
                sel_y=2'b01;
                sel_z=2'b01;
                sel_T=2'b01;
                en_T=1'b1;
                en_z=1'b1;
                done=1'b0;
                nxt_state=S2;
                	// Write your output and next state equations here
                end
            S2 :
                  begin
				sel_x=2'b10;
                sel_y=2'b10;
                sel_z=2'b10;
                sel_T=2'b10;
                en_T=1'b1;
                en_z=1'b1;
                done=1'b0;
                nxt_state=S3;
                	// Write your output and next state equations here
                end
              S3 :
                begin
				sel_x=2'b11;
                sel_y=2'b11;
                sel_z=2'b10;
                sel_T=2'b11;
                en_T=1'b1;
                en_z=1'b0;
                done=1'b0;
                nxt_state=S4;
                	// Write your output and next state equations here
                end
			// Define the rest of the states
                S4 :
                begin
				
                sel_x=2'b11;
                sel_y=2'b11;
                sel_z=2'b11;
                sel_T=2'b11;
                en_T=1'b0;
                en_z=1'b1;
                done=1'b0;
                nxt_state=S5;
                
                	// Write your output and next state equations here
                end

                S5 :
                begin
				sel_x=2'b00;
                sel_y=2'b00;
                sel_z=2'b00;
                sel_T=2'b00;
                en_T=1'b1;
                en_z=1'b1;
                
                done=1'b1;
                
                	// Write your output and next state equations here
                end
            default: 
                begin

                sel_x=2'b00;
                sel_y=2'b00;
                sel_z=2'b00;
                sel_T=2'b00;
                en_T=0;
                en_z=0;
                done=1'b0;
                nxt_state=S0;    
				// Don't forget the default
                end            
        endcase
        
    end

endmodule


module reg_with_enable #(parameter N = 32) (clk, rst, en, X, O );
    input [N:0] X;
    input clk;
    input rst;
    input en;
    output [N:0] O;
    
    reg [N:0] R;
    
    always@(posedge clk) begin
        if (rst) begin
            R <= {N{1'b0}};
        end
        if (en) begin
            R <= X;
        end
    end
    assign O = R;
endmodule







/*-------------------Supporting Modules--------------------*/
/*------------- Iterative Karatsuba: 32-bit Karatsuba using a single 16-bit Module*/

module mult_16(X, Y, Z);
input [15:0] X;
input [15:0] Y;
output [31:0] Z;

assign Z = X*Y;

endmodule




module full_adder(a, b, cin, S, cout);
input a;
input b;
input cin;
output S;
output cout;

assign S = a ^ b ^ cin;
assign cout = (a&b) ^ (b&cin) ^ (a&cin);

endmodule


module check_subtract (A, B, C);
 input [7:0] A;
 input [7:0] B;
 output [8:0] C;
 
 assign C = A - B; 
endmodule



/* N-bit RCA adder (Unsigned) */
module adder_Nbit #(parameter N = 32) (a, b, cin, S, cout);
input [N-1:0] a;
input [N-1:0] b;
input cin;
output [N-1:0] S;
output cout;

wire [N:0] cr;  

assign cr[0] = cin;


generate
    genvar i;
    for (i = 0; i < N; i = i + 1) begin
        full_adder addi (.a(a[i]), .b(b[i]), .cin(cr[i]), .S(S[i]), .cout(cr[i+1]));
    end
endgenerate    


assign cout = cr[N];

endmodule


module Not_Nbit #(parameter N = 32) (a,c);
input [N-1:0] a;
output [N-1:0] c;

generate
genvar i;
for (i = 0; i < N; i = i+1) begin
    assign c[i] = ~a[i];
end
endgenerate 

endmodule


/* 2's Complement (N-bit) */
module Complement2_Nbit #(parameter N = 32) (a, c, cout_comp);

input [N-1:0] a;
output [N-1:0] c;
output cout_comp;

wire [N-1:0] b;
wire ccomp;

Not_Nbit #(.N(N)) compl(.a(a),.c(b));
adder_Nbit #(.N(N)) addc(.a(b), .b({ {N-1{1'b0}} ,1'b1 }), .cin(1'b0), .S(c), .cout(ccomp));

assign cout_comp = ccomp;

endmodule


/* N-bit Subtract (Unsigned) */
module subtract_Nbit #(parameter N = 32) (a, b, cin, S, ov, cout_sub);

input [N-1:0] a;
input [N-1:0] b;
input cin;
output [N-1:0] S;
output ov;
output cout_sub;

wire [N-1:0] minusb;
wire cout;
wire ccomp;

Complement2_Nbit #(.N(N)) compl(.a(b),.c(minusb), .cout_comp(ccomp));
adder_Nbit #(.N(N)) addc(.a(a), .b(minusb), .cin(1'b0), .S(S), .cout(cout));

assign ov = (~(a[N-1] ^ minusb[N-1])) & (a[N-1] ^ S[N-1]);
assign cout_sub = cout | ccomp;

endmodule



/* n-bit Left-shift */

module Left_barrel_Nbit #(parameter N = 32)(a, n, c);

input [N-1:0] a;
input [$clog2(N)-1:0] n;
output [N-1:0] c;


generate
genvar i;
for (i = 0; i < $clog2(N); i = i + 1 ) begin: stage
    localparam integer t = 2**i;
    wire [N-1:0] si;
    if (i == 0) 
    begin 
        assign si = n[i]? {a[N-t:0], {t{1'b0}}} : a;
    end    
    else begin 
        assign si = n[i]? {stage[i-1].si[N-t:0], {t{1'b0}}} : stage[i-1].si;
    end
end
endgenerate

assign c = stage[$clog2(N)-1].si;

endmodule



