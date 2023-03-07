typedef enum
{
/* https://handwiki.org/wiki/List_of_instruction_sets */

/* 1. Advanced Digital Chips Inc. */
EISC, // (AE32K) 32-bit embedded core

/* 2. Altera (Later, Intel) */
Nios, // II:[2] 32-bit RISC, CPU core optimized for implementation in Altera FPGAs

/* 3. AMD */
AM29000, // (112 Instructions):[3] RISC
/* AMD extensions to x86
    3DNow! (21 instructions):[4] An extension for floating-point arithmetic
    3DNow! Extensions (5 instructions):[5]:316 An extension for 3DNow!.
    ABM - Advanced Bit Manipulation
    SSE4a - streaming SIMD extensions 4a
    SSE5 - streaming SIMD extensions 5 (170 instructions, proposal)
    F16C - FP16 conversion operations, a revision of part of the proposed SSE5
    XOP - eXtended Operations, a revision of part of the proposed SSE5
    TBM - Trailing Bit Manipulation
*/
AMD64, // :[5][6][7][8][9][10] 64-bit extension of x86, originally named x86-64, adopted by Intel

/* 4. Analog Devices, Inc. (ADI) */
SHARC, // [11]
Blackfin, // [12]

/* 5. Apollo Computer Inc. */
Apollo_PRISM, // : 32-bit VLIW RISC (Not to be confused with DEC Prism or m88k.)

6. ARC International (Later, Synopsys)
ARC:[13] Pre-configurable core
ARCompact
ARCtangent

7. Arm
ARM architecture[14]
ARMv1
ARMv2
ARMv3
ARMv4
ARMv5
Extensions
Thumb[15]
DSP[16]
Jazelle
VFPv2 - vector floating point
ARMv6
Extensions
Thumb-2[17]
TrustZone[18]
SIMD[19]
ARMv7
Extensions
Thumb-2
NEON - media acceleration technology[20]
VFPv3
ARMv8-A
Extensions

8. AT&T (Later, Lucent Then Agere Then LSI, Then Avago and Intel)
For StarCore DSP architecture, refer to Motorola section.

DSP16
DSP16A[21]
DSP1600
ATT2100[22]
WE32K (WE 32000, WE 32100, WE 32200)[23]

9. Atmel
ATiny
AVR[24]
AVR32[25]

10. Axis Communications
ETRAX CRIS

11. Burroughs (Later, Unisys)
B1000 series
B2000/B3000/B4000 series
B5000 series and B6000/B7000 series

12. Cambridge Consultants / Cambridge Silicon Radio (Later, CSR Plc)
XAP series

13. CDC (Control Data Corporation)
CDC 160 series
CDC 924
CDC 1604
CDC 3000 24-bit
CDC 3000 48-bit
CDC 6000 series
CDC 7600
CDC STAR-100
Cyber 70 and 170 series
CDC Cyber 180
CDC Cyber 200

14. Commonwealth Scientific and Industrial Research Organisation
CSIRAC

15. Cray Research, Inc. (Later, Silicon Graphics, Inc., Then Cray Inc.)
Cray-1 supercomputer

16. Cypress
M8C Core: 8-bit MCU PSoC 1

17. Data General
Nova:[26] 16-bit CISC
ECLIPSE:[27] 16-bit CISC
ECLIPSE/MV:[28] 32-bit CISC

18. DEC (Digital Equipment Corporation)
PDP-1: 18-bit CISC minicomputer
PDP-4/PDP-7/PDP-9/PDP-15: 18-bit CISC minicomputer
PDP-5/PDP-8/PDP-12: 12-bit CISC minicomputer
PDP-6/PDP-10/DECSYSTEM-20: 36-bit CISC mainframe
PDP-11:[29] 16-bit CISC minicomputer
VAX:[30] 32-bit CISC
Prism: 32-bit RISC
Alpha: 64-bit RISC

19. Donald Knuth
Introduced in the textbook of Prof. Donald Knuth

MIX
MMIX[31]

20. DSP Group (Later, CEVA, Inc.)
DSP Group and Parthus Technologies plc were merged into CEVA, Inc. in 2002.

Oak DSP Core
Teak Series
Teak DSP Core
TeakLite DSP Core
CEVA-TeakLite-4
CEVA-X
CEVA-XC
CEVA-XC4000
CEVA-XM4

21. Eckert-Mauchly Computer Corporation (Later, Remington Rand Then Sperry Then Unisys)
UNIVAC 1: The first commercial computer produced in the United States

22. Elliott Brothers
Elliott Automation

23. EnSilica
eSi-RISC

24. Fairchild
Clipper[32]

25. Fujitsu (Later, Cypress)
FR Series: 32-bit RISC MCU
FR30/FR60 Family ()[33]
FR80/FR81 Family ():[34] Some instructions are removed and added against FR30 and FR60.
FR-V: VLIW and vector processor based RISC
F2MC Series
F2MC-16 Family ():[35] 16-bit MCU
F2MC-8L/F2MC-8FX Family ():[36] 8-bit MCU

26. General Electric (Later, Honeywell, Then Honeywell/Bull, and Then NEC Corporation)
GE-200 series:[37] Small main frame, 20-bit word machine
GE-400 series:[38] Middle mainframe, 24-bit word machine
GE-412:[39] 20-bit word machine
GE-600 series/Honeywell 6000 series: Large main frame, 36-bit CISC, word machine, LSB on left
GE-625/635[40]
GE-645:[41] Multics available
Toshiba TOSBAC-5600: GECOS-3 and ACOS-6 available
HIS (Honeywell Information Systems) 6025, 6030, 6040, 6050, 6060, 6070, 6080: GCOS available
HIS 6180: Multics available
HIS Series 60 Level 66 and Level 66/DPS: GCOS available
HIS Series 60 Level 68 and Level 68/DPS: Multics available[42]
HIS DPS-8: GCOS available
HIS DPS-8M: Multics available[42]
Honeywell Bull DPS-88: GCOS available
Honeywell Bull DPS-8000:[43] GCOS available
NEC ACOS Series 77 System 600, 700, 600S, 800, 900:[44] Succeeded Toshiba's business. GCOS-3 and ACOS-6 available
NEC ACOS System 1000 (HIS DPS90), 2000 (HIS DPS/9000), 3900 (HIS DPS/9000-900):[45] ACOS-6 and GCOS-8 available, No Multics

27. General Instrument Microelectornics (Later, Microchip Technology Incorporated)
The company was established as a subsidiary of General Instrument in 1987, then became an independent company as Microchip Technology in 1989.

CP1600: 16-bit microprocessor
SP0256 - Speech processor[46]
PIC microcontroller
Mid-range PIC[47]
PIC16[48]
PIC1605: NMOS 8-bit microcontroller, the basis of PIC architecture
PIC17[48]
PIC18[48]
dsPIC30F[49]
dsPIC33
28. Hennessy (,Prof.) and Patterson (,Prof.)
DLX:[50] Introduced as educational-use ISA in their famous textbooks; "Computer architecture : a quantitative approach" and "Computer organization and design : the hardware/software interface." GNU assembler is available.
Stanford MIPS: Basis of MIPS architecture by Prof. John L. Hennessy
Berkeley RISC: Basis of SPARC architecture by Prof. David Patterson
29. Hewlett-Packard
HP 2100
FOCUS
HP 3000 "Classic" CISC
PA-RISC
PA-RISC 1.0
PA-RISC 1.1
MAX-1 SIMD extensions
PA-RISC 2.0
MAX-2 SIMD extensions
30. Hitachi (Later, Renesas)
HD6309 ():[51] An extension for Motorola MC6809
HD64180: Z80-based embedded MCU
H8 Family
H8/500 (63 instructions)[52]
For other H8 Series, refer to Renesas section.
SuperH RISC engine Family[53]
SuperH, Inc. SH-5 ():[54][55] 64-bit RISC
For other SH Series, refer to Renesas section.
31. Holtek Semiconductor
HT RISC:[56]:24-36 8-bit RISC MCU
32. Honeywell
These are instruction sets introduced by Honeywell; for the instruction sets from General Electric, refer to the General Electric section.

Datamatic 1000,[57] H-400, H-1400, H-800,[58] H-1800,[59] and H-1800-II:[60] 48-bit word machine with 3 address format
Series 200 model 200/1200/2200:[61] A character-oriented two-address commercial computer
Honeywell Model 8200:[62] A system containing a word-processing subsystem based on the H-800 and a Variable Length Field (VLF) processor based on the H-200
DDP Series 16 model 316 and 516:[63] 16-bit minicomputer
33. IBM
IBM 1130/IBM 1800
IBM 1400 series/IBM 7010
IBM 1620/IBM 1710
IBM 37xx
IBM 3790
IBM 650
IBM 701
IBM 704/IBM 709/IBM 7090/IBM 7094/IBM 7040/IBM 7044
IBM 702/IBM 705/IBM 7080
IBM 7070/IBM 7072/IBM 7074
IBM 7030 Stretch
System/360 (32-bit CISC) and successors
System/370: 32-bit CISC
System/390: 32-bit CISC
z/Architecture: 64-bit CISC
IBM 8100
IBM Series/1
IBM System/3[64]
IBM System/4 Pi
AP-101: Used in the moon flights
IBM System/32[65]
IBM System/34[66]
IBM System/36[67]
IBM System/38/IBM AS/400/IBM System i MI code[68][69]
IBM System/7
ROMP
Power Architecture
POWER ISA: POWER1, the RISC Single Chip, POWER2
PowerPC ISA: POWER3
PowerPC v2.00: POWER4
PowerPC v2.01:[70] POWER5
PowerPC v2.02:[71] POWER5+
Book E (Enhanced PowerPC) by NXP,[72] by STM[73]
Cell Broadband Engine Architecture (Including the PowerXCell 8i):[74]:51 PowerPC ISA v2.02 + 8× SPE (Synergistic Processing Elements) for vector/SIMD multimedia extensions
PowerPC AS
Power ISA
Power ISA v2.03:[75] POWER6
Power ISA v2.04
Power ISA v2.05: POWER6+
Power ISA v2.06B:[76] POWER7
Power ISA v2.07B (for POWER8 & POWER8 with Nvidia NVLink)[77][78]
Power ISA v3.0B (for POWER9)[77][79]
PowerQUICC by NXP: PowerPC + plural of QUICC vector processor elements
34. Infineon Technologies AG
C166[80]
C500[81]
TriCore[82][83]
35. INMOS and XMOS
by Prof. David May

Transputer[84]
XCore[85]
36. Intel
4004 (46 instructions)[86]:288/621
4040 (60 instructions)[86]:218/621
8008 (48 instructions)[86]:420/621
8080 (111 Instructions),[87] 8085 (113 Instructions)[88]
8021 (66 Instructions)[89]
8022 (73 Instructions)[90]
MCS-41 (also known as 8041) (87 instructions)[91]
MCS-48 (also known as 8048) (93 instructions)[92]
MCS-51 (also known as 8051) (111 instructions)[93]:§3[94]:§2[95]
iAPX 432[96][97][98][99]
i860:[100][101] 32/64-bit VLIW RISC
i960 (also known as 80960) (FIX MI core instructions with 11 addressing modes):[102][103] 32-bit RISC
IA-64 (also known as Itanium):[104] Originated at Hewlett-Packard (HP), and later jointly developed by HP and Intel
x86, See: x86 instruction listings
8086/8088, 80186/80188, 80286: 16-bit CISC
IA-32:[105][106] 32-bit CISC
x86-64: 64-bit extension of x86, originally developed by AMD as AMD64[6][7][8][9][10]
Intel 64, Intel's version of x86-64
Extensions[107]
FPU (x87) - Floating-point-unit (FPU) instructions
MMX - MMX SIMD instructions
MMX Extended - extended MMX SIMD instructions
3DNow! (21 instructions):[4] An extension for floating-point arithmetic (AMD)
3DNow! Extensions (5 instructions):[5]:316 An extension for 3DNow! (AMD)
SSE - streaming SIMD extensions (SSE) instructions (70 instructions)
SSE2 - streaming SIMD extensions 2 instructions (144 new instructions)
SSE3 - streaming SIMD extensions 3 instructions (13 new instructions)
SSSE3 - supplemental streaming SIMD extensions (16 instructions)
SSE4.1 - streaming SIMD extensions 4, Penryn subset (47 instructions)
SSE4.2 - streaming SIMD extensions 4, Nehalem subset (7 instructions)
SSE4 - All streaming SIMD extensions 4 instructions (both SSE4.1 and SSE4.2)
SSE4a - streaming SIMD extensions 4a (AMD)
SSE5 - streaming SIMD extensions 5 (170 instructions, proposal from AMD)
F16C - FP16 conversion operations (from AMD), a revision of part of the proposed SSE5
XOP - eXtended Operations (AMD), a revision of part of the proposed SSE5
ABM - Advanced Bit Manipulation (from AMD)
TBM - Trailing Bit Manipulation (AMD)
XSAVE - XSAVE instructions
AVX - advanced vector extensions instructions
FMA - fused multiply-add instructions
AES - Advanced Encryption Standard instructions
CLMUL - Carry-less mtiply (PCLMULQDQ) instruction
Cyrix - Cyrix-specific instructions
AMD - AMD-specific instructions (older than K6)
SMM - System management mode instructions
SVM - Secure virtual machine instructions
PadLock - VIA PadLock instructions
37. IPFlex (Later, Tokyo Keiki, Inc.)
DAP/DNA-2: A Dynamic Reconfigurable Processor, jointly developed by IPFlex and Fujitsu.
38. Lattice Semiconductor
LatticeMico8[108]
LatticeMico32[109]
39. Lebedev Institute of Precision Mechanics and Computer Engineering
BESM-1
BESM-2
BESM-3, BESM-4
BESM-6
40. Maxim Integrated
MAXQ[110]:§18
41. MIPS Technologies
MIPS architecture
MIPS I
MIPS II
MIPS III
MIPS IV[111]
MIPS V
MIPS16
MIPS32
MIPS64
MDMX
Loongson Technology
Loongson is a Chinese company. In its earlier stage, its architecture was MIPS like because of patents problem until a deal in 2007. In 2011, it formally licensed MIPS32 and MIPS64.
Loongson 1: 32-bit MIPS like. Lacking 4 instructions by patent issue.
Loongson 3: MIPS64 quad core. Over 200 instructions are added for x86 emulation.
42. MIT's Lisp Machine (Later, Symbolics, Inc., Then Symbolics)
CADR Lisp Machine
43. Mitsubishi Electric (Later, Renesas)
Mitsubishi D10V[112]
Mitsubishi D30V[113]
740 - 8-bit 6502 superset
7700 Family - 16-bit
7700 Series (103 instructions)
7751 Series (109 instructions)
7900 Series (203 instructions):[114] 16-bit MCU
44. MOS Technology (Later, Commodore Semiconductor Group)
MOS/CSG 6502:[115] 8-bit CISC
MOS/CSG 65CE02: added extra registers and instructions, having 6502 emulation
Western Design Center 65816: 16-bit CISC, having 6502 emulation
45. Motorola (Later, Freescale and Then NXP Semiconductor)
For the Power Architecture, refer to IBM section.

6800 Family: 8-bit CISC
Motorola 6800 (107 instructions)[116]
Motorola 6801 (98 instructions)[117]
Motorola 6805 (86 instructions)[118]
Motorola 6809 (94 instructions)[119]
Hitachi HD6309: Having some extensions.
Freescale HC11 (62 instructions)
CPU16 Family: 16-bit CISC
Freescale HC16
68000 Family: 32-bit CISC
68000[120]
68010[121]
68020 and 68030
68881 and 68882 FPUs
68040
68060
683XX
ColdFire
88000: 32-bit RISC
DSP56800[122]
StarCore DSP Architecture:[123]:23 Jointly developed by Agere and Motorola.
SC100(V1)/SC140(V2):[124] Baseline
SC140e(V3): Minor additions
SC3400(V5): Video, SIMD2
SC3850(V6): Control ISA, Dual MPY, Cache instructions
SC3900FP:[125] SIMD8
46. National Semiconductor
COP8
CR16
NSC800[126]
NS320xx[127][128][129]
47. NCR Corporation
48. NEC Corporation
Supercomputer
SX architecture: A Scalar Processing Unit + eight Vector Processing Elements
SX-1 and SX-2
SX-3
SX-4 - SX-9[130][131][132][133]
Earth Simulator
SX-ACE[134]
Mainframe:
Large mainframe (refer to General Electric section)
Middle mainframe
NEC ACOS Series 77 System 300, 400, 500:[44] ACOS-4 OS, 32-bit byte machine
NEC ACOS System 1500[135]
Small mainframe
NEC ACOS Series 77 System 200:[44] ACOS-2 OS, byte machine
Minicomputer
MS-4100: 32-bit CISC, byte machine, big endian
49. NEC's Semiconductor Operations (Later, NEC Electronics, Then Renesas Electronics)
μCOM-1600, also known as μPD768 (93 instructions): 16-bit single-chip CPU released in 1978. Used for NEC System 50 office computer.[136][137]
17K Famiy ():[138]:§3:575/690-665/690 4-bit MCU
75 Family
μPD7500: 4-bit MCU
μPD7500 set ? (106 instructions)
μPD7500 Set A (92 instructions)[139]:52-58
μPD7500 Set ? (67 instructions)
75X Series ():[140]:§10[141]:§10 4-bit MCU
75XL Series (): 4-bit MCU
μCOM-87 Family
μCOM-87 Series (): 8-bit MCU
μCOM-87LC Series (): 8-bit MCU
μCOM-87AD Series
NMOS version (158 instructions):[142]:§13 8-bit MCU
CMOS version (159 instructions):[143]:§14 8-bit MCU
78K Family
78K/0 Series: 8-bit MCU. Refer to Renesas section
78K0R Series: 16/8-bit MCU. Refer Renesas RL78 section
78K/0S Series: 8-bit MCU. Refer to Renesas section
78K/1 Series(64 instructions):[144]:§18 8-bit MCU
78K/2 Series(65 instructions):[145]:§18[146]:§20 8-bit MCU
78K/3 Series(111 instructions with macro service):[147]:§18 16/8-bit MCU
78K/4 Series(113 instructions with macro service):[148] 16/8-bit MCU
78K/6 Series ( with macro service): 16-bit MCU
V60/V70, V80 (119, 123 instructions):[149] 32-bit CISC, little endian
V810/V830:[150] 32-bit RISC, little endian
V850: 32-bit RISC, Refer to Renesas section
50. OpenRISC Community
OpenRISC 1000:[151] ORBIS32; 32-bit RISC, big endian
51. Parallax, Inc.
Propeller P8X32A:[152] 32-bit RISC
52. PEZY Computing
PEZY-1
PEZY-SC
53. RCA
CDP1802[153]
Spectra 70 (System/360 compatible in user mode ("problem state"), not compatible in kernel mode ("supervisor state"))
54. Renesas
The semiconductor operations of Hitachi and Mitsubishi Electric were transferred to Renesas Technology Corporation on April 1, 2003. In addition, NEC Electronics Corporation, a subsidiary of NEC Corporation, and Renesas Technology were merged into Renesas Electronics Corporation on April 1, 2010.

RL78 Family:[154] 8/16-bit CISC MCU, similar ISA to 8-bit 78K/0 legacy CISC, accumulator-based architecture, 2 operand instructions, 1-5 byte non-uniform length instructions, 13 addressing modes, non-orthogonal instruction set, little endian, 3-stage pipeline
RL78-S1 Core (74 instructions):[154] 8-bit ALU, 8× 8-bit registers, No register-banks
RL78-S2 Core (75 instructions):[154] 16-bit ALU, 8× 8-bit registers, 4× register-banks
RL78-S3 Core (81 instructions):[154] 16-bit ALU, multiplication/division/MAC extension
78K0R (80 instructions divided into 15 groups):[155] 75 instructions out of 80 are identical to RL78.
16-bit ALU, 8× 8-bit registers, 4× register-banks,
RX Family:[156]:14PDF8 32-bit CISC MCU, general-purpose-register-based architecture, 2 and 3 operand instructions, 1-8 byte non-uniform length instructions, 16× 32-bit GPRs, highly orthogonal instruction set, bi-endian, optional single precision floating-point arithmetic operations, 5-stage synchronous pipeline
RXv1 (90 instructions/10 addressing modes, basic:73, float:8, DSP:9)[157]
RXv2 (109 instructions/11 addressing modes, basic:75, float:11, DSP:23)[158]
RH850 Family: Upward compatible with V850 ISA, disclosed to Automotive industry customers only
V850 Family: 32-bit RISC MCU, general-purpose-register-based architecture, load/store architecture, 2 operand instructions, basically 2-byte and 4-byte 2-way form instructions (having extension), relatively orthogonal instruction sets, branch with interlock, 32× 32-bit GPRs, little endian, optional single and double precision floating-point arithmetic operations, 5- or 7-stage synchronous pipeline
V850 (74 instructions):[159] Simple but sined-load/store, saturation arithmetic
V850E (81 instructions)[160] Unsigned-load/store extension, CISCy extension
V850E1 (80 (83) instructions)[161]
V850ES (80 instructions)[162]
V850E1F (96 instructions): Floating-point arithmetic extension for single precision
(Also known as PHOENIX-Fx)
V850E2 (89 instructions):[163] Superscaler
V850E2S (98 instructions):[164] Memory protection
V850E2M (98+α instructions):[165] Memory protection, floating-point for single and double precision
V850E3 (): SIMD extension, loop extension with branch predictor
SuperH RISC engine Family:[53] 32/64-bit RISC MCU/MPU, general-purpose-register-based architecture, load/store architecture, 2-byte uniform length instruction set, relatively orthogonal instruction sets, branch with delay slots, 16× 32-bit GPRs with partially 2 banks of 8 registers, 1× 32-bit global base register, 2× 32-bit MAC register, 1× 32-bit procedure register, optional 2 banks of 16× 32-bit floating-point registers, optional 2× 40-bit and 6× 32-bit DSP registers, bi-endian, 5- or 7-stage synchronous pipeline
SH-1 (56 instructions):[166] 32-bit RISC
SH-2 (62 instructions):[166] 32-bit RISC
SH2-DSP (154 instructions):[166] 32-bit RISC
SH-2A (91 instructions):[167] 32-bit RISC
SH2A-FPU (112 instructions):[167] 32-bit RISC
SH-3 (68 instructions):[168] 32-bit RISC
SH-3E (84 instructions):[168] 32-bit RISC
SH3-DSP (160 instructions):[168] 32-bit RISC
SH-4 (93 instructions):[169] 32-bit RISC
SH4A (103 instructions):[170] 32-bit RISC
SH4AL-DSP (226 instructions):[171] 32 -bit RISC
For SH-5 Series, refer to Hitachi section
78K Family:
78K0 Series (48 instructions):[172] 8-bit MCU, accumulator-based architecture, 8× 8-bit registers, 4× register-banks, non-pipelined
78K0S Series (47 instructions):[173] 8-bit simplified version of 78K0, no mul/div insn., no register-bank, etc.[174]
For 78K/1, 78K/2, 78K/3, 78K/4, and 78K/6 Series, refer to NEC's semiconductor section.
For 78K0R, refer to RL78 section.
R8C Family: 16-bit CISC MCU
R8C/Tiny Series (89 instructions)[175]
M16C Family (106 instructions):[176] 16-bit CISC
M32C Family (108 instructions):[177] 16/32-bit CISC MCU, 2 banks of 4 × 16-bit data and 2 × 24-bit address registers
M32R Family (83 instructions):[178] 32-bit RISC MCU
M32R-FPU (100 instructions):[179] floating-point arithmetic extension
H8SX Family ():[180] 32-bit MCU
H8S Family ():[181] 32-bit MCU
H8 Family
H8/300 Series ():[182] 16/8-bit MCU
H8/300H Series ():[183] 16/8-bit MCU
H8/300L Series ():[184] 16/8-bit MCU
For H8/500 Series, refer to Hitachi section
720 Family (135 instruction):[185]:83-139 4-bit MCU, accumulator
740 Family (71 instructions):[186] A 8-bit 6502 superset
55. Rockwell Collins
56. Samsung Electronics
SAM8
57. Scenix Semiconductor (Later, Ubicom Then Qualcomm)
Ubicom was acquired by Qualcomm in 2012.

SX Series: 8-bit microcontroller
IP Series: 32-bit microprocessor (IP2000, IP8000, etc.)
58. Signetics
2650 (75 instructions)[187]
59. SpaceWire UK
Raptor-16:[188] 16-bit CISC
60. STMicroelectronics (Formerly, SGS-Thomson)
For SPC5 Power Architecture Book E product line, refer to IBM section.

STM8 (80 instructions, 20 addressing modes):[189] 8-bit MCU
ST10 ("FIX-ME" basic instructions and "FIX-ME" MAC instructions):[190] 16-bit MCU
ST40 ():[191][192] 32-bit RISC SuperH family SH-4 architecture, jointly developed with Hitachi
61. Sun Microsystems (Later, Oracle)
SPARC[193][194]
SPARC-V7[195]
SPARC-V8[195]
SPARC-V9:[195] Ultra-SPARC I (1995)
Extensions[195]:97
VIS 1: UltraSPARC I (1995)
VIS 2: UltraSPARC III (2001)
EMAf: SPARC T3 (2010)
VIS 3, VIS 3B: SPARC T3 (2010)
IMA: SPARC T4 (2011)
SPARC4 (Core S3, 2011)
Java bytecode for Java VM:[196] 32-bit stack machine, RISC + CISC, big endian
picoJava:[197] Having extended Bytecodes
JEM1[198] by Rockwell Collins: A super-set of JVM Bytecode
JEMCore-II by aJile Systems, Inc. (aJ-102 and aJ-200[199]:107): Having another extended Bytecodes
Java Card VM Bytecode:[200] 16-bit stack machine, having completely different ISA from Java VM
MAJC: VLIW
62. Tensilica (Later, Cadence)
CPUs
Xtensa:[201] Pre-configurable architecture
DSPs
HiFi Audio and Voice DSPs
Vision DSPs
ConnX D2[202]
ConnX BBE16[202]
63. Texas Instruments
9900[203]
9940[204]
9980[205]
MSP430[206]
TMS320 series
64. UNIVAC (Later, Unisys)
1100/2200 series
65. University of California, Berkeley (UCB) (Later, RISC-V Foundation)
RISC-V:[207] Open ISA introduced in 2010
66. University of Cambridge
EDSAC: The first practical stored-program computer
CAP computer
67. University of Texas at Austin, and University of Illinois at Urbana-Champaign
LC-3: 16-bit RISC ISA for educational use
LC-3b:[208] Modified variant introduced with hardware microarchitecture
68. University of Tokyo
TAC:[209] A tube computer developed in 1959
69. U.S. Military
ENIAC: One of the earliest electronic general-purpose computer
MIL-STD-1750A:[210] The U.S.'s military standard computer, 16-bit RISC
Apollo Guidance Computer: Used in the moon flights
70. Xerox
71. Xilinx
MicroBlaze[211]
PicoBlaze[212]
72. Zilog (Later, a Subsidiary of IXYS Corporation)
Z80, Z800, Z280, Z180, Z380, eZ80[213]
Z8, eZ8[214]
Z8000, Z80000[215]
} Emulator_Arch_type ;

typedef struct
{
    const 
} Emulator_Arch ;

typedef struct
{
    const  Type
} Emulator_Opcode_Struct ;

#if Emulator_RISCV
const Emulator_Opcode_Struct 
#endif /* Emulator_RISCV */