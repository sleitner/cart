erase
lweight 3
ltype 0
ctype 0


location 5000 32000 4000 32000


ticksize 0 0 0 0
notation -1 2 -1 2


define IBOX 1


macro plot 2 {
	da $1
	read { lr1 1 a1 2 }

	set r = 10**lr1 if ( a1 > 0 )
	set a = a1 if ( a1 > 0 )

	if ( $IBOX ) {
		limits 0 6 0 1.5
		expand 1.5
		box
		expand 2.0
		xlabel r
		ylabel q_{F}
		define IBOX 0
	}

	if ( $2 ) {
		ptype 30 0
		set q = a/as
		points r q
	} else {
		set as = a
	} 

}

plot OUT/prof_te.res 0
plot OUT/prof_ss.res 1


macro fit 1 {
	set r2 = r**2
	set s1 = r2/(r2+$1)
	set w = ( r2 < 1 ) ? 1.306 : s1
	connect r w
}


