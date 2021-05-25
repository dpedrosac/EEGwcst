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
c1z	 1	C	1	4	1639	   ;	# 16. COL
e1a	 2	C	2	3	1605	   ;	
c4v	 3	C	3	2	1725	   ;	
c4r	 4	C	4	1	1903	   ;	
b2a	 5	C	5	3	1541	   ;		
e3r	 6	N	1	3	1787	   ;	# 17. NUM
c4v	 7	N	2	4	1785	   ;	
c2z	 8	N	3	2	1757	   ;	
e4r	 9	N	4	4	1909	   ;	
e4a	10	N	5	4	1926	   ;	
t2a	11	F	1	1	1943	   ;	# 18. FOR
c1v	12	F	2	3	1845	   ;	
e3r	13	F	3	2	1565	   ;	
b2r	14	F	4	4	1557	   ;	
t3z	15	F	5	1	1685	   ;		

};
