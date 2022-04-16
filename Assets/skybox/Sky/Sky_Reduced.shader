Shader "Custom/SkyReduced"
{
		//planetary shadow ref. https://www.youtube.com/watch?v=CnXFe4MVkkU
	Properties
	{
		_UpDir("Up Direction", Vector) = (0.0, 1.0, 0.0)	
		_DirProfile("RGB direction profile for particles",2D) = "white" {}
		_AdensZero("density",float)=1.0
		_altitude("altitude",float)=1.0
		_sunAngularSize("sun rad angular",Range(0,360))=1.0
	//	_scatteringCol("Scattering Color",Color ) = (1,1,1,1)
	//	_atmosCol("AtmosColor",Color ) = (1,1,1,1)
	//	_reflectCol("ReflectColor",Color ) = (1,1,1,1)
	//	_absorbRatio("absorb to reflect ratio",Range(0,0.5))=0.1
		_plantetRad("Planetary Radius KM",float)=6000.0
		_atmosH("Atmosphere Height KM",float)=11.0

	}
	SubShader
	{
		Tags
		{
			"PreviewType" = "Plane"
		}
		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#include "UnityCG.cginc"

			struct VertIn
			{
				float4 vertex : POSITION;
			};
	
			float3 _UpDir;
			float4 _LightColor0;
			sampler2D _DirProfile;
			float _AdensZero;
			float _altitude;
			float _sunAngularSize;
			float _absorbRatio;
			float4 _atmosCol;
			float4 _reflectCol;
			float4 _scatteringCol;
			float _plantetRad;
			float _atmosH;

			struct VertOut
			{
				float4 position : SV_POSITION;
				float4 worldPosition : TEXCOORD0;
			};

			VertOut vert(VertIn v)
			{
				VertOut o;

				o.position = UnityObjectToClipPos(v.vertex);
				o.worldPosition = mul(unity_ObjectToWorld, v.vertex);

				return o;
			}


			float4 frag(VertOut i) : SV_Target
			{
				//cam to pixel ray vector.
				float3 fragDir = normalize(i.worldPosition - _WorldSpaceCameraPos);
				float3 upDir = normalize(_UpDir);
				float elevation = dot(fragDir, upDir); // -1..1
				float atmosDensityMask = 1-clamp(elevation,0,1); //1 at horizon 0 at top
				elevation = (1 + elevation) / 2; // 0..1 bottom..top
				
				float3 dirToLight = normalize(_WorldSpaceLightPos0);

				//	float4 skyColor = SkyColor(camLightHMatch, elevation);

				float lightFacingNormalized = dot(fragDir, dirToLight); // -1..1
				lightFacingNormalized = (1 + lightFacingNormalized) / 2; // 0..1

				//float4 lightColor = LightColor(lightFacingNormalized);
				//base formula for angular capture is Arccos(sqrt(d^2-R^2)/d)
				//we also must divide by 2 to get radius
				/*
				TODO: move to sky controller. dont need these calculations in shader.
				float sunVisAngle = 1.0-(acos(sqrt(_AdensZero*_AdensZero-_altitude*_altitude)/_AdensZero)/3.1415926538)/2.0;
   				float testangle=(180.0/360);
				float sunangle=(0.5/360)/2;
				*/
				float sunMask = 1-(step(lightFacingNormalized,1-(_sunAngularSize/360)/2));
				float4 sunColor=_LightColor0;
				float4 sun = sunColor*sunMask;
				
				//Planetary shadow mask
				float horizonDist=200;
				float3 domePoint = horizonDist*fragDir;
				
				float3 lightHzDir = normalize(float3(dirToLight.x,0,dirToLight.z));//normalize(dirToLight*float3(1,0,1));
				float3 shadowPlaneOffset=float3(0,-50,0);
				float3 lightHzPoint = lightHzDir*horizonDist+shadowPlaneOffset; //shifting plane below visual horizon.

				float planetShadow =1;
				if(dirToLight.y<0){
					float3 planarVec=float3 (-lightHzDir.z,0,lightHzDir.x); //this works.
					float3 shadowPlaneUp = (cross(dirToLight,planarVec));
					float planarDistance = -dot(shadowPlaneUp,(domePoint-lightHzPoint));
					planetShadow=(planarDistance/horizonDist)*8;
					planetShadow=clamp(planetShadow,0,1);
				}

				float Rp = _plantetRad; //radius of planet km. earth = 6000
				float Ha = _atmosH; //atmosphere max height km. earth = 12/ venus=200

				float sk=0.99;
				float3 Lpower=4000000*sunColor.rgb;
				float eyeMaxLight=4000;
				float3 EarthAD0=float3(0.007,0.02,0.02)*2; //AD is percentage of particle interactions over 1KM.
				float3 FogAD0=0.02;
				float3 Air=float3(0,0,0.013);
				float3 Dust=float3(0.003,0.008,0);
				
				float3 Adens0 = _AdensZero; //real density
				float3 Adens0SK = 1.3;// kg/m3 default density at zero level which is used for sk
// sk - scatter koefficient
// prk - propagation coefficient.

				float3 skRGB=float3(0.0012, 0.0019, 0.0032)*13; //how much light is scattered by 1000m3 at adens0
				float3 prkRGB=float3(0.6, 0.4, 0.2); //percent of scattered light that propagates along with light ray. in integrals it is used as sk*prk. it is percent of scatter.
				//prkRGB=float3(0.8, 0.8, 0.8);

				///prkRGB=0.9;

				float3 focusRGB;//amonut of light that scatters to the eye from point on viewline
				float viewang=max(0.5*(dot(fragDir,dirToLight)+1),0);
				focusRGB=tex2D(_DirProfile,float2(1-viewang,0.2));

				float viewAlt=_altitude; // km
				float Lmax=0; //distance to atmos ray exit from inside
				float LmaxOuter=0; //distance to atmos ray entrance from outside.
				float3 viewPoint=float3(0,Rp+viewAlt,0); //in coordinates where center of planet is at 0,0,0
				float rmax = Rp+Ha;
				float rmax2=rmax*rmax;

				{
					float dtoProj=dot((float3(0,0,0)-viewPoint),fragDir); //good!
					float3 projPt = viewPoint+dtoProj*fragDir; // works!
				//	if (viewAlt>Ha)
				//		Lmax = dtoProj-sqrt(rmax2-pow(length(projPt),2)); //km
				//		else
						Lmax = sqrt(rmax2-pow(length(projPt),2))+dtoProj; //km 
						LmaxOuter = dtoProj-sqrt(rmax2-pow(length(projPt),2)); //km 

						Lmax=max(Lmax,0);
						LmaxOuter=max(LmaxOuter,0);
				}

				float surfaceMask=0;
				{
					//km, as planet radius
					float arg=Rp/(Rp+viewAlt);
					float ylim = sqrt(1-arg*arg);// or cos(asin(arg));
					//for negative altitude make -ylim
					surfaceMask= step(fragDir.y+ylim,0);
				}

				float surfdist=0;

				float LmaxEye=Lmax-LmaxOuter;

				if(surfaceMask>0){
					float dtoProj=dot((float3(0,0,0)-viewPoint),fragDir); //good!
					float3 projPt = viewPoint+dtoProj*fragDir; // works!
					surfdist = dtoProj-sqrt(Rp*Rp-pow(length(projPt),2));//-dtoProj; //km
					LmaxEye=surfdist-LmaxOuter;
					Lmax=surfdist;
				}
				// this Lmax could be easily used as fast approximation of atmosphere gradient. for cheap effects.
				//return LmaxEye/1000;



				// probe of depth using ray. page: Atmos.RAY)
				float atmosDepth=0;
				float3 endPoint=viewPoint+Lmax*fragDir;
				float integralDepth=0;
				{
					float rV=length(viewPoint);
					float rE=length(endPoint);
					//basically extending area by lmax/ha)
					//accurate
					//integralDepth=(((rE-rV)*(rE-2*Rp+rV-2*Adens0*Ha))/(2*Ha))*(Lmax/(Ha));
					integralDepth=(Adens0*Ha/2)*(Lmax/Ha); //as simple as it gets.
					//float h0=0; //height at point 0; and height at endpoint is always 1;
					//integralDepth=((Adens0*h0*h0)/2*Ha)+((Adens0*Ha/2)-Adens0*h0)*(Lmax/Ha);
					//formula for reflection.
					float maxdens=10;
					float angleCond=(1-((1-clamp(dot(fragDir,dirToLight),0.8,1))/0.2));
					float reflected=angleCond*(integralDepth/50);
					if (reflected>1)reflected =0;
					float absorbed=integralDepth/50;
					//return reflected+max(sun-absorbed*2,0);
				}
				
				// lensing effect on sun.
				float refrSun=lerp(_sunAngularSize,_sunAngularSize*1.2,integralDepth/150);
				sunMask=1-(step(lightFacingNormalized,1-(refrSun/360)/2));
	


				float dust=0;
				//dust reflection  https://www.youtube.com/watch?v=eDqsY68hNBs
				// https://www.youtube.com/watch?v=4AwNFxPFQBM

				float rayScatToEye=0;
				{
					float r = Rp+Ha;
					float r2=r*r;
					float3 scanstart=viewPoint+LmaxOuter*fragDir;
					float3 scanend=scanstart+LmaxEye*fragDir;//viewPoint+fragDir*LmaxOuter;
					
					float distToProjs=dot((float3(0,0,0)-scanstart),dirToLight); //good! distance from scanpoint to projection on light vector
					float distToProje=dot((float3(0,0,0)-scanend),dirToLight); //good!

					float3 projPts = scanstart+distToProjs*dirToLight; 
					float3 projPte = scanend+distToProje*dirToLight; // works!
					
					float endlPv= sqrt(r2-pow(length(projPts),2))+distToProjs;
					float endlPve = sqrt(r2-pow(length(projPte),2))+distToProje; //distance from top viewpoint to light exit
					endlPv=max(endlPv,0);
					endlPve=max(endlPve,0);

					float3 endpPv = scanstart+dirToLight*endlPv;
					float3 endpPve = scanend+dirToLight*endlPve;



					float RayH=dot(fragDir*Lmax,float3(0,1,0)); //and this is good density approximation;
					//(1-endlPve/(Rp*2)) is a variant for atmospheric shadow due to ray absorption.

					float hAdjBot=length((scanstart+endpPv)*0.5)-Rp; 
					float hAdjTop=length((scanend+endpPve)*0.5)-Rp;

					float3 roTop=Adens0-(hAdjTop/Ha)*Adens0; //average density for ray at the bottom
					float3 roBot=Adens0-(hAdjBot/Ha)*Adens0;
					float3 roEye=Adens0/2;
	
					//working formula without fade from point to eye. and improper density
					//float3 rsRGB=Adens0*Lpower*Lmax*focusRGB*skRGB*(8* pow(((prkRGB-1)*roTop*skRGB+1),endlPve)+7*pow(((prkRGB-1)*roBot*skRGB+1),endlPv))/30;
					//proper one with powers

					float3 toEyeSK=skRGB/1000; // SK IS USED FOR 1000M
					//skRGB/1000 was;
					float3 rsRGB=Adens0*Lpower*Lmax*focusRGB*(toEyeSK)*(8* pow(((prkRGB-1)*skRGB+1),endlPve*roTop/Adens0SK)+7*pow(((prkRGB-1)*skRGB+1),endlPv*roBot/Adens0SK))/30;

					
					float3 eyeFade=pow((1+skRGB*(prkRGB-1)),(LmaxEye*roEye/Adens0SK)/2);
						rsRGB=rsRGB*eyeFade; //eyefade is due to distance from scattering particle to eye.
					
					float4 remSun= float4(sunMask*Lpower*pow((1+skRGB*(prkRGB-1)),Lmax*roBot/Adens0SK)/eyeMaxLight,1);

					return ( float4(rsRGB,1)/eyeMaxLight + remSun*(1-surfaceMask));//*(1-surfaceMask);
				}
			}
			ENDCG
		}
	}
}