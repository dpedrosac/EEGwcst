# WCST TRADICIONAL

scenario = "WIS1";
scenario_type = trials;

active_buttons = 4;
button_codes = 31, 32, 33, 34;   
target_button_codes = 31, 32, 33, 34;
write_codes = true;
pulse_width = 50 ; 
screen_width = 800; 
screen_height = 600; 
screen_bit_depth = 16; 


default_trial_type = first_response;      
default_font_size = 25;
default_font = "Arial black";       
response_logging = log_active ;    
# default_all_responses = false;

default_background_color = 255, 255, 255; # white
default_text_color = 0, 0, 0; # black

# define variables
$w = 177 ; 		# initial width of bipmap files (original 177)
$h = 177 ; 		# initial height of bipmap files (original 177)

$t1= 1000 ;    
$t2=  700 ;
$st = 1 ;      # Neuroscan stimulus code          
$c = "N" ;



begin;
box { height = 200; width = 200; color = 255,255,255; } wblanck;
bitmap { filename = "keys.bmp"; width = 393; height = 97; } wkeys;
bitmap { filename = "b1a.bmp"; width = $w; height = $h; } wb1a;
bitmap { filename = "b1v.bmp"; width = $w; height = $h; } wb1v;
bitmap { filename = "b2a.bmp"; width = $w; height = $h; } wb2a;
bitmap { filename = "b2r.bmp"; width = $w; height = $h; } wb2r;
bitmap { filename = "b3r.bmp"; width = $w; height = $h; } wb3r;
bitmap { filename = "b3v.bmp"; width = $w; height = $h; } wb3v;
bitmap { filename = "c1v.bmp"; width = $w; height = $h; } wc1v;
bitmap { filename = "c1z.bmp"; width = $w; height = $h; } wc1z;
bitmap { filename = "c2r.bmp"; width = $w; height = $h; } wc2r;
bitmap { filename = "c2z.bmp"; width = $w; height = $h; } wc2z;
bitmap { filename = "c4r.bmp"; width = $w; height = $h; } wc4r;
bitmap { filename = "c4v.bmp"; width = $w; height = $h; } wc4v;
bitmap { filename = "e1a.bmp"; width = $w; height = $h; } we1a;
bitmap { filename = "e1z.bmp"; width = $w; height = $h; } we1z;
bitmap { filename = "e3r.bmp"; width = $w; height = $h; } we3r;
bitmap { filename = "e3z.bmp"; width = $w; height = $h; } we3z;
bitmap { filename = "e4a.bmp"; width = $w; height = $h; } we4a;
bitmap { filename = "e4r.bmp"; width = $w; height = $h; } we4r;
bitmap { filename = "t2a.bmp"; width = $w; height = $h; } wt2a;
bitmap { filename = "t2z.bmp"; width = $w; height = $h; } wt2z;
bitmap { filename = "t3v.bmp"; width = $w; height = $h; } wt3v;
bitmap { filename = "t3z.bmp"; width = $w; height = $h; } wt3z;
bitmap { filename = "t4a.bmp"; width = $w; height = $h; } wt4a;
bitmap { filename = "t4v.bmp"; width = $w; height = $h; } wt4v;

picture {} default;
picture { box wblanck; x = 0; y = 0; } blanck;
picture { bitmap wkeys; x = 0; y = 130; bitmap wb1a; x = 0; y = 0; } b1a;
picture { bitmap wkeys; x = 0; y = 130; bitmap wb1v; x = 0; y = 0; } b1v;
picture { bitmap wkeys; x = 0; y = 130; bitmap wb2a; x = 0; y = 0; } b2a;
picture { bitmap wkeys; x = 0; y = 130; bitmap wb2r; x = 0; y = 0; } b2r;
picture { bitmap wkeys; x = 0; y = 130; bitmap wb3r; x = 0; y = 0; } b3r;
picture { bitmap wkeys; x = 0; y = 130; bitmap wb3v; x = 0; y = 0; } b3v;
picture { bitmap wkeys; x = 0; y = 130; bitmap wc1v; x = 0; y = 0; } c1v;
picture { bitmap wkeys; x = 0; y = 130; bitmap wc1z; x = 0; y = 0; } c1z;
picture { bitmap wkeys; x = 0; y = 130; bitmap wc2r; x = 0; y = 0; } c2r;
picture { bitmap wkeys; x = 0; y = 130; bitmap wc2z; x = 0; y = 0; } c2z;
picture { bitmap wkeys; x = 0; y = 130; bitmap wc4r; x = 0; y = 0; } c4r;
picture { bitmap wkeys; x = 0; y = 130; bitmap wc4v; x = 0; y = 0; } c4v;
picture { bitmap wkeys; x = 0; y = 130; bitmap we1a; x = 0; y = 0; } e1a;
picture { bitmap wkeys; x = 0; y = 130; bitmap we1z; x = 0; y = 0; } e1z;
picture { bitmap wkeys; x = 0; y = 130; bitmap we3r; x = 0; y = 0; } e3r;
picture { bitmap wkeys; x = 0; y = 130; bitmap we3z; x = 0; y = 0; } e3z;
picture { bitmap wkeys; x = 0; y = 130; bitmap we4a; x = 0; y = 0; } e4a;
picture { bitmap wkeys; x = 0; y = 130; bitmap we4r; x = 0; y = 0; } e4r;
picture { bitmap wkeys; x = 0; y = 130; bitmap wt2a; x = 0; y = 0; } t2a;
picture { bitmap wkeys; x = 0; y = 130; bitmap wt2z; x = 0; y = 0; } t2z;
picture { bitmap wkeys; x = 0; y = 130; bitmap wt3v; x = 0; y = 0; } t3v;
picture { bitmap wkeys; x = 0; y = 130; bitmap wt3z; x = 0; y = 0; } t3z;
picture { bitmap wkeys; x = 0; y = 130; bitmap wt4a; x = 0; y = 0; } t4a;
picture { bitmap wkeys; x = 0; y = 130; bitmap wt4v; x = 0; y = 0; } t4v;


sound { wavefile { filename = "t1000Hz.wav"; } ; } good;
sound { wavefile { filename = "t0500Hz.wav"; } ; } bad;  


trial {
	trial_duration = forever ;
	picture {
		text { caption = 
  "Instructions 1. 
  >>> Press a key to start practice <<< ";
		};
		x = 0; y = 0;
	};
	time = 0;
};                     

trial {                         
	sound bad;
	duration = 250 ;
	code = "fb_shi3" ;     
	port_code = 10 ;
	time = $t1 ;
} shift3;

trial {      
   sound bad;
	duration = 250 ;
	code = "fb_shi2" ;  
	port_code = 20 ;
	time = $t1 ;
} shift2;

trial {                         
	sound bad;
	duration = 250 ;
	code = "fb_fail" ;
	port_code = 50 ;
	time = $t2 ;
} failed;
 
trial {      
   sound good;
	duration = 250 ;
	code = "fb_ant" ;       
	port_code = 40 ;
	time = $t2 ;
} anti;

trial {      
   sound good;
	duration = 250 ;
	code = "fb_rep1" ;    
	port_code = 21 ;
	time = $t2 ;
} repeat1;  

trial {      
   sound good;
	duration = 250 ;
	code = "fb_rep2" ;    
	port_code = 22 ;
	time = $t1 ;
} repeat2;  

trial {      
   sound good;
	duration = 250 ;
	code = "fb_rep3" ;    
	port_code = 23 ;
	time = $t1 ;
} repeat3;  

trial {      
   sound good;
	duration = 250 ;
	code = "fb_rep4" ;    
	port_code = 24 ;
	time = $t2 ;
} repeat4;  

trial {      
   sound good;
	duration = 250 ;
	code = "fb_rep5" ;    
	port_code = 25 ;
	time = $t2 ;
} repeat5;  



TEMPLATE "mcst.tem" {
stim	tr	c	st	tb	t1	      ;	# WCST TRADICIONAL
e3r	1	N	1	3	1711	   ;	# 1. NUM
t4v	2	N	2	4	1847	   ;	
t2z	3	N	3	2	1948	   ;	
b3r	4	N	4	3	1948	   ;	
e1a	5	N	5	1	1529	   ;	
c1z	6	N	6	1	1524	   ;	
t2a	7	N	7	2	1938	   ;	
e3z	8	N	7	3	1780	   ;	
c4r	9	N	7	4	1943	   ;	
b2a	10	F	1	4	1865	   ;	# 2. FOR
t2z	11	F	2	1	1637	   ;	
c1v	12	F	3	3	1738	   ;	
e4r	13	F	4	2	1538	   ;	
t4a	14	F	5	1	1604	   ;	
c2z	15	F	6	3	1719	   ;	
b3v	16	F	7	4	1981	   ;	
c2r	17	F	7	3	1776	   ;	
c1v	18	C	1	2	1851	   ;	# 3. COL
t2a	19	C	2	3	1529	   ;	
t3z	20	C	3	4	1726	   ;	
e4a	21	C	4	3	1582	   ;	
c4r	22	C	5	1	1886	   ;	
t4a	23	C	6	3	1650	   ;	
b2r	24	C	7	1	1923	   ;	
t3v	25	F	1	1	1926	   ;	# 4. FOR
b3r	26	F	2	4	1554	   ;	
c2z	27	F	3	3	1613	   ;	
t4a	28	F	4	1	1651	   ;	
c4v	29	F	5	3	1530	   ;	
t3z	30	F	6	1	1882	   ;	
e3z	31	F	7	2	1775	   ;	
b1v	32	F	7	4	1624	   ;	
b2a	33	N	1	2	1698	   ;	# 5. NUM
e3z	34	N	2	3	1941	   ;	
c4r	35	N	3	4	1903	   ;	
e1z	36	N	4	1	1723	   ;	
t3z	37	N	5	3	1557	   ;	
c2r	38	N	6	2	1825	   ;	
t3v	39	N	7	3	1527	   ;	
c1z	40	C	1	4	1548	   ;	# 6. COL
t3v	41	C	2	2	1938	   ;	
b1a	42	C	3	3	1629	   ;	
e3r	43	C	4	1	1774	   ;	
e1z	44	C	5	4	1502	   ;	
b3r	45	C	6	1	1628	   ;	
t2z	46	C	7	4	1647	   ;	
b3v	47	N	1	3	1588	   ;	# 7. NUM
b1a	48	N	2	1	1626	   ;	
c1v	49	N	3	1	1871	   ;	
t4v	50	N	4	4	1653	   ;	
c2z	51	N	5	2	1988	   ;	
e4r	52	N	6	4	1643	   ;	
t4a	53	N	7	4	1860	   ;	
b3r	54	N	7	3	1587	   ;	
c1z	55	N	7	1	1662	   ;	
b1v	56	C	1	2	1899	   ;	# 8. COL
b3v	57	C	2	2	1531	   ;	
c1z	58	C	3	4	1886	   ;	
b2a	59	C	4	3	1539	   ;	
e1a	60	C	5	3	1716	   ;	
e4r	61	C	6	1	1868	   ;	
t4v	62	C	7	2	1796	   ;	
c2r	63	F	1	3	1659	   ;	# 9. FOR
c4r	64	F	2	3	1922	   ;	
t4v	65	F	3	1	1801	   ;	
b1a	66	F	4	4	1648	   ;	
e1z	67	F	5	2	1776	   ;	
e1a	68	F	6	2	1566	   ;	
b1v	69	F	7	4	1506	   ;	
t2z	70	C	1	4	1551	   ;	# 10. COL
e1z	71	C	2	4	1766	   ;	
c2r	72	C	3	1	1686	   ;	
b1a	73	C	4	3	1919	   ;	
c4v	74	C	5	2	1524	   ;	
b3v	75	C	6	2	1815	   ;	
b3r	76	C	7	1	1665	   ;	
t4a	77	C	7	3	1612	   ;	
t3z	78	N	1	3	1699	   ;	# 11. NUM
t2a	79	N	2	2	1657	   ;	
e3z	80	N	3	3	1769	   ;	
b2a	81	N	4	2	1747	   ;	
t3v	82	N	5	3	1750	   ;	
b2r	83	N	6	2	1957	   ;	
b3v	84	N	7	3	1866	   ;	
t4a	85	N	7	4	1992	   ;	
e1a	86	N	7	1	1714	   ;	
c4v	87	F	1	3	1788	   ;	# 12. FOR
t3z	88	F	2	1	1629	   ;	
e3z	89	F	3	2	1736	   ;	
t4v	90	F	4	1	1988	   ;	
b1a	91	F 	5	4	1824		;
e4a	92	F	6	2	1689	   ;	
b2r	93	F	7	4	1564	   ;	
t2z	94	F	7	1	1922	   ;	
c1v	95	C	1	2	1970	   ;	# 13. COL
t3v	96	C	2	2	1557	   ;	
t3z	97	C	3	4	1794	   ;	
e4r	98	C	4	1	1761	   ;	
b1v	99	C	5	2	1560	   ;	
e3z	100	C	6	4	1615	   ;	
b2r	101	C	7	1	1739	   ;	
b1v	102	N	1	1	1957	   ;	# 14. NUM
e4a	103	N	2	4	1973	   ;	
e1z	104	N	3	1	1504	   ;	
c4v	105	N	4	4	1783	   ;	
b2r	106	N	5	2	1511	   ;	
b1a	107	N	6	1	1659	   ;	
c2z	108	N	7	2	1769	   ;	
b1v	109	F	1	4	1872	   ;	# 15. FOR
e4r	110	F	2	2	1565	   ;	
b3v	111	F	3	4	1689	   ;	
e4a	112	F	4	2	1628	   ;	
c4r	113	F	5	3	1503	   ;	
b2a	114	F	6	4	1718	   ;	
t3v	115	F	7	1	1767	   ;	
c1z	116	F	7	3	1672	   ;	
c1z	117	C	1	4	1639	   ;	# 16. COL
e1a	118	C	2	3	1605	   ;	
c4v	119	C	3	2	1725	   ;	
c4r	120	C	4	1	1903	   ;	
b2a	121	C	5	3	1541	   ;	
t4v	122	C	6	2	1996	   ;	
e4a	123	C	7	3	1706	   ;	
e3r	124	N	1	3	1787	   ;	# 17. NUM
c4v	125	N	2	4	1785	   ;	
c2z	126	N	3	2	1757	   ;	
e4r	127	N	4	4	1909	   ;	
e4a	128	N	5	4	1926	   ;	
t3z	129	N	6	3	1717	   ;	
c1v	130	N	7	1	1573	   ;	
t2a	131	F	1	1	1943	   ;	# 18. FOR
c1v	132	F	2	3	1845	   ;	
e3r	133	F	3	2	1565	   ;	
b2r	134	F	4	4	1557	   ;	
t3z	135	F	5	1	1685	   ;	
c1z	136	F	6	3	1645	   ;	
e3r	137	F	7	2	1951	   ;	

};
