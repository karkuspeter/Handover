package application;

import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import static com.kuka.roboticsAPI.motionModel.BasicMotions.ptp;

import com.kuka.common.ThreadUtil;
import com.kuka.roboticsAPI.applicationModel.RoboticsAPIApplication;
import com.kuka.roboticsAPI.controllerModel.Controller;
import com.kuka.roboticsAPI.deviceModel.LBR;
import com.kuka.roboticsAPI.geometricModel.CartDOF;
import com.kuka.roboticsAPI.geometricModel.Frame;
import com.kuka.roboticsAPI.geometricModel.ObjectFrame;
import com.kuka.roboticsAPI.geometricModel.PhysicalObject;
import com.kuka.roboticsAPI.geometricModel.Tool;
import com.kuka.roboticsAPI.geometricModel.Workpiece;
import com.kuka.roboticsAPI.geometricModel.World;
import com.kuka.roboticsAPI.motionModel.ISmartServoRuntime;
import com.kuka.roboticsAPI.motionModel.LIN;
import com.kuka.roboticsAPI.motionModel.ServoMotion;
import com.kuka.roboticsAPI.motionModel.SmartServo;
import com.kuka.roboticsAPI.motionModel.controlModeModel.CartesianImpedanceControlMode;
import com.kuka.roboticsAPI.sensorModel.ForceSensorData;
import com.kuka.roboticsAPI.uiModel.ApplicationDialogType;
import com.kuka.roboticsAPI.userInterface.ServoMotionUtilities;
import com.kuka.generated.ioAccess.*;


public class OptitrackHandoverAndras_test extends RoboticsAPIApplication {
	private Controller kuka_Sunrise_Cabinet_1;
	private LBR lbr_iiwa_14_R820_1;
	private Tool ToolGripper;
	private Tool WaterBottle;

	private Workpiece Object_05;

	private CartesianImpedanceControlMode mode = null ; 
	private CartesianImpedanceControlMode modeHard = null ; 
	private CartesianImpedanceControlMode modeHandShake = null ; 

	private GripperIOGroup FingerTip;
	private ATIIOGroup ForceTorqueSensor ;

	private UdpClient client = null;

	Boolean _IsRunning = true;
	double[] OptiTrackData;
	int count = 0;

	double robotBaseX = 0.0;
	double robotBaseY = 0.0;
	double robotBaseZ = 0.0;

	int projectNumber = 99;
	int bodyOrHand = 100;
	int fingerControlOn = 100;
	int handoverDialogChoice =99;
	
	int CommandedFingerPosition = 255;
	int CommandedFingerPositionPrev = 150;
	double HumanForce= 0.0;
	int timeframe ; 
	double RB1_x = 0.0 ;
	double RB1_y = 0.0 ;
	double RB1_z = 0.0 ;
	double RB1_qx = 0.0 ;
	double RB1_qy = 0.0 ;
	double RB1_qz = 0.0 ;
	double RB1_qw = 0.0 ;
	
	double posXmsr;
	double posYmsr;
	double posZmsr ;


	double RB2_x = 0.0 ;
	double RB2_y = 0.0 ;
	double RB2_z = 0.0 ;
	double RB2_qx = 0.0 ;
	double RB2_qy = 0.0 ;
	double RB2_qz = 0.0 ;
	double RB2_qw = 0.0 ;
	
	double qx = 0.0;
	double qy = 0.0;
	double qz = 0.0;
	double qw = 0.0;

	double bodyDistance = 10000.0;
	double handDistance = 10000.0;
	double bodyDistancePrev = 10000.0;
	double bodyVelocity = 0.0;
	double dt = 1.0/120.0; // 120 FPS
	double handTrackDistanceLimit = 200.0;
	
	boolean followHand = false;

	long currTime = (long) 0.0;
	long prevTime = (long)0.0;
	long dtLoop = (long)0.0;
	
	double rotAlpha = 0.0;
	double rotBeta = 0.0;
	double rotGamma = 0.0;

	
	@Override
	public void initialize()
	{
		// Locate the "first" Lightweight Robot in the system
		lbr_iiwa_14_R820_1 = ServoMotionUtilities.locateLBR(getContext());
		// FIXME: Set proper Weights or use the plugin feature

		kuka_Sunrise_Cabinet_1 = (Controller) getContext().getControllers().toArray()[0];

		FingerTip= new GripperIOGroup(kuka_Sunrise_Cabinet_1);
		//	ForceTorqueSensor = new ATIIOGroup(kuka_Sunrise_Cabinet_1);

		ToolGripper=getApplicationData().createFromTemplate("ToolGripper");		
		WaterBottle = getApplicationData().createFromTemplate("WaterBottle");

		mode = new CartesianImpedanceControlMode() ; 
		mode.parametrize(CartDOF.A).setStiffness(5) ;  	// rotation about Z
		mode.parametrize(CartDOF.B).setStiffness(5);		// rotation about Y
		mode.parametrize(CartDOF.C).setStiffness(5);		// rotation about X

		mode.parametrize(CartDOF.X).setStiffness(300) ;
		mode.parametrize(CartDOF.Y).setStiffness(300) ;
		mode.parametrize(CartDOF.Z).setStiffness(300) ;
		mode.parametrize(CartDOF.ALL).setDamping(0.7) ;
	}


	public int ReceiveData()
	{

		String recvStr = client.getString();
		String[] optiTrackStr = recvStr.split(" ");
		OptiTrackData=new double[optiTrackStr.length];
		for(int j=0;j<optiTrackStr.length;j++)
		{
			OptiTrackData[j]=Double.valueOf(optiTrackStr[j]);

		}
		return 1;
	}

	@Override 
	public void dispose(){

		_IsRunning = false ; 
		System.out.println(" closing sockets in Dispose Block"); 

		if(client != null)
		{
			client.kill();
		}

		super.dispose();
	}

	/**
	 * Move to an initial Position WARNING: MAKE SHURE, THAT the pose is
	 * collision free.
	 */
	private void moveToInitialPosition()
	{
		lbr_iiwa_14_R820_1.move(ptp( Math.toRadians(0.0), Math.toRadians(66.0), Math.toRadians(10.0), Math.toRadians(-82.0) , Math.toRadians(13.0), Math.toRadians(-66.0), Math.toRadians(40.0)).setJointVelocityRel(0.25));

		if (FingerTip.getActReq() != 9)
		{
			FingerTip.setActReq(9);
		}
		FingerTip.setSpeed(255);
		FingerTip.setForce(100);

		FingerTip.setPosReq(255); // close the gripper

		ThreadUtil.milliSleep(2000); // wait for 2 seconds
	}

	private double[] getForce() {

		ForceSensorData torData = lbr_iiwa_14_R820_1.getExternalForceTorque(ToolGripper.getFrame("/WaterBottleCP"));

		double ForceX = torData.getForce().getX();
		double ForceY = torData.getForce().getY();
		double ForceZ = torData.getForce().getZ();

		double ForceMagnitude = Math.sqrt(	Math.pow(ForceX, 2) + 
				Math.pow(ForceY, 2) + 
				Math.pow(ForceZ, 2)		);

		double[] Forces = new double[3];
		Forces[0] = ForceX;
		Forces[1] = ForceY;
		Forces[2] = ForceZ;
		Forces[3] = ForceMagnitude;

		return Forces;
	}

	public void run() {

		try{
			try {
				client = new UdpClient();
				client.start();

				ReceiveData();

				if(OptiTrackData[0]==-1)
				{
					System.err.println("No optitrack data");
					return;
				}
				else
				{
					String S_reply = "Connected";
					System.out.println(S_reply);
				}

			} 

			catch (Exception ex) {
				System.err.println("Can't accept client connection. ");
			}	

			ToolGripper.attachTo(lbr_iiwa_14_R820_1.getFlange());

			moveToInitialPosition();

			
			
			
			if (projectNumber > 5){
				projectNumber = HandoverDialog();
			}
			
			if (projectNumber == 1){ // Hand following
				bodyOrHand = TrackingDialog();
			}
			
			if (projectNumber == 4){
				while (fingerControlOn != 2) {
					fingerControlOn = 	FingerControlDialog();
					if (fingerControlOn == 0) {//opening
						FingerTip.setPosReq(0);
					} else if (fingerControlOn == 1) { // closing 
						FingerTip.setPosReq(255);
					} else if (fingerControlOn == 2) {
						WaterBottle.attachTo(ToolGripper.getFrame("/WaterBottleCP"));
						break;
					}
				}
			}
			
			if (projectNumber == 0) { // Handover
				
				while (handoverDialogChoice != 3){
					handoverDialogChoice = DialogInitHandover();
					if (handoverDialogChoice == 0) {//opening
						FingerTip.setPosReq(0);
					} else if (handoverDialogChoice == 1) { // closing and attaching
						FingerTip.setPosReq(255);
					} else if (handoverDialogChoice == 2) {
						ToolGripper.getFrame("/WaterBottleCP").move(ptp(getApplicationData().getFrame("/PreObjectPosition")).setJointVelocityRel(.4));
						FingerTip.setPosReq(0);
						ThreadUtil.milliSleep(1500);
						ToolGripper.getFrame("/WaterBottleCP").move(ptp(getApplicationData().getFrame("/ObjectPosition")).setJointVelocityRel(.4));
						FingerTip.setPosReq(255);
						WaterBottle.attachTo(ToolGripper.getFrame("/WaterBottleCP"));
						ThreadUtil.milliSleep(1500);
						ToolGripper.getFrame("/WaterBottleCP").move(ptp(getApplicationData().getFrame("/PreObjectPosition")).setJointVelocityRel(.4));
						ToolGripper.getFrame("/WaterBottleCP").move(ptp(getApplicationData().getFrame("/HandoverPosition")).setJointVelocityRel(.4));
						break;
					} else if (handoverDialogChoice == 3){
						WaterBottle.attachTo(ToolGripper.getFrame("/WaterBottleCP"));
						break;
					}
				}
			}
			
			SmartServo aSmartServoMotion = new SmartServo(lbr_iiwa_14_R820_1.getCurrentJointPosition());
			aSmartServoMotion.useTrace(true);

			// for Automatic mode 0.25, for T1 mode 1
			aSmartServoMotion.setJointAccelerationRel(0.4);
			aSmartServoMotion.setJointVelocityRel(0.2);
			
			aSmartServoMotion.setMinimumTrajectoryExecutionTime(5e-3);

			if ( projectNumber == 0){ // Handover with object
				if (aSmartServoMotion.validateForImpedanceMode(WaterBottle) != true)
				{
					getLogger().info("Validation for SmartServo Compliant control failed");
				}

				WaterBottle.getRootFrame().moveAsync(aSmartServoMotion.setMode(mode));
			} else { // no object
				if (aSmartServoMotion.validateForImpedanceMode(ToolGripper) != true)
				{
					getLogger().info("Validation for SmartServo Compliant control failed");
				}

				ToolGripper.getFrame("/WaterBottleCP") .moveAsync(aSmartServoMotion.setMode(mode));
			}

			// Fetch the Runtime of the Motion part
			ISmartServoRuntime theServoRuntime = aSmartServoMotion.getRuntime();

			while(_IsRunning)
			{

				ReceiveData();

				if(OptiTrackData[0]==-1)
				{
					System.out.println(" stop motion ! at " + theServoRuntime.toString());
					// Stop the motion
					theServoRuntime.stopMotion();

					break;
				}
				// Synchronize with the realtime system
				theServoRuntime.updateWithRealtimeSystem();


				if (projectNumber == 3){ // check optitrack data
					count+=1;

					RB1_x = (OptiTrackData[2] - robotBaseX)*1000;
					RB1_y = (OptiTrackData[3] - robotBaseY)*1000;
					RB1_z = (OptiTrackData[4] - robotBaseZ)*1000;

					System.out.println("RB1 XYZ = " + String.valueOf( RB1_x ) + ", "+ String.valueOf( RB1_y ) + ", "+ String.valueOf( RB1_z ));

					RB2_x = (OptiTrackData[10] - robotBaseX)*1000;
					RB2_y = (OptiTrackData[11] - robotBaseY)*1000;
					RB2_z = (OptiTrackData[12] - robotBaseZ)*1000;

					System.out.println("RB2 XYZ = " + String.valueOf( RB2_x ) + ", "+ String.valueOf( RB2_y ) + ", "+ String.valueOf( RB2_z ));
					System.out.println("Entered the loop " + String.valueOf(count) + " times");
					break;
				}else if (projectNumber == 2) // context estimation
				{

					count +=1;

					RB2_x = (OptiTrackData[10] - robotBaseX)*1000;
					RB2_y = (OptiTrackData[11] - robotBaseY)*1000;
					RB2_z = (OptiTrackData[12] - robotBaseZ)*1000;

					currTime = System.nanoTime();

					bodyDistance = Math.sqrt( Math.pow( (RB2_x), 2) 
							+ 	Math.pow( ( RB2_y), 2) 
							+ 	Math.pow( ( RB2_z), 2) 	)  ;

					bodyVelocity = Math.abs(bodyVelocity * 0.5  + 0.5 * (bodyDistance-bodyDistancePrev)/dt);
					dtLoop =  (currTime-prevTime);

					bodyDistancePrev = bodyDistance;
					if (bodyDistance <= 2000.0) {

						dtLoop = (currTime-prevTime)/(long)1e9;

						getLogger().info("Estimated body velocity is: " + String.valueOf(bodyVelocity) + ", @ " + String.valueOf(RB2_x) + "/" + String.valueOf(RB2_y) + "/" + String.valueOf(RB2_z) + ", distance= " + String.valueOf(bodyDistance) + ", dt = " + String.valueOf(dt)+ ", dtLoop = " + String.valueOf(dtLoop) + ", entered the loop: " + String.valueOf(count) + " times");

						FingerTip.setPosReq(150);
						break;
					}
					prevTime= currTime;


				} else if (projectNumber == 1) //hand following no object
				{
					count +=1;

					RB1_x = (OptiTrackData[2] - robotBaseX)*1000;
					RB1_y = (OptiTrackData[3] - robotBaseY)*1000;
					RB1_z = (OptiTrackData[4] - robotBaseZ)*1000;

					RB2_x = (OptiTrackData[10] - robotBaseX)*1000;
					RB2_y = (OptiTrackData[11] - robotBaseY)*1000;
					RB2_z = (OptiTrackData[12] - robotBaseZ)*1000;

					currTime = System.nanoTime();

					bodyDistance = Math.sqrt(  		Math.pow( (RB2_x), 2) 
							+ 	Math.pow( ( RB2_y), 2) 
							+ 	Math.pow( ( RB2_z), 2) 	)  ;
					
					handDistance = Math.sqrt(  		Math.pow( (RB1_x), 2) 
							+ 	Math.pow( ( RB1_y), 2) 
							+ 	Math.pow( ( RB1_z), 2) 	)  ;

					if (bodyDistance <= 1500.0) {
						followHand = true;
					} else {followHand = false;}

					if (RB1_z > 1000.0) {
						getLogger().info("Hand raised => finish hand tracking");
						theServoRuntime.stopMotion();
						break;
					}

					if (followHand) {FingerTip.setPosReq(50); 

					Frame goalFrame = new Frame();  
					goalFrame = lbr_iiwa_14_R820_1.getCurrentCartesianPosition(ToolGripper.getFrame("/WaterBottleCP"));

					if ( Math.abs(RB1_z) >= 800.0) {RB1_z = 800.0;}
					if ( Math.abs(RB1_z) <= -100.0) {RB1_z = -100.0;}

					if (bodyOrHand == 1){
						// compute the absolute motion w.r.t. World 
						goalFrame.setX(RB1_x);
						goalFrame.setY(RB1_y);
						goalFrame.setZ(RB1_z);

						double theta1 = Math.atan(RB1_y/RB1_x);
						goalFrame.setAlphaRad(77.0*Math.PI/180.0 + theta1);
						goalFrame.setBetaRad(0);
						goalFrame.setGammaRad(90*Math.PI/180.0);
					
						} else {
							goalFrame.setX(RB2_x);
							goalFrame.setY(RB2_y);
							goalFrame.setZ(RB2_z);

							double theta1 = Math.atan(RB2_y/RB2_x);
							goalFrame.setAlphaRad(77.0*Math.PI/180.0 + theta1);
							goalFrame.setBetaRad(0);
							goalFrame.setGammaRad(90*Math.PI/180.0);
							
							
						}

					theServoRuntime.setDestination(goalFrame, World.Current.getRootFrame());
					

					}
					else {FingerTip.setPosReq(150); 
					}

				
				} else if (projectNumber == 5) { // Orientation test
					
					RB1_x = (OptiTrackData[2] - robotBaseX)*1000;
					RB1_y = (OptiTrackData[3] - robotBaseY)*1000;
					RB1_z = (OptiTrackData[4] - robotBaseZ)*1000;

					qx = OptiTrackData[5];
		  			qy = OptiTrackData[6];
		  			qz = OptiTrackData[7];
		  			qw = OptiTrackData[8];
		  			
		  			rotAlpha = Math.atan2(-2*(qx*qy-qw*qz), qw*qw +qx*qx -qy*qy-qz*qz ) *180/Math.PI;
		  			rotBeta = Math.asin(2*(qx*qz + qw*qy))*180/Math.PI;
		  			rotGamma = Math.atan2(-2*(qy*qz-qw*qx), qw*qw - qx*qx -qy*qy +qz*qz)*180/Math.PI;
				
		  			getLogger().info("Pos: " + String.valueOf(RB1_x) + ", " + String.valueOf(RB1_y) + ", " + String.valueOf(RB1_z) +
		  					"Quat: " + String.valueOf(qx) + ", " + String.valueOf(qy) + ", " + String.valueOf(qz) + ", " + String.valueOf(qw) +
		  					"Angles: " + String.valueOf(rotAlpha) + ", " + String.valueOf(rotBeta) + ", " + String.valueOf(rotGamma));
		  			ThreadUtil.milliSleep(500);
				} else if (projectNumber == 4) { //fixed point handover
// Checking force data					
//					ForceSensorData torData = lbr_iiwa_14_R820_1.getExternalForceTorque(
//							ToolGripper.getFrame("/WaterBottleCP"),
//							ToolGripper.getFrame("/WaterBottleCP"));
//
//					double ForceX = torData.getForce().getX() ;
//					double ForceY = torData.getForce().getY() ;
//					double ForceZ = torData.getForce().getZ();
//					
//					getLogger().info("Forces: " + String.valueOf(ForceX) + ", " + String.valueOf(ForceY) + ", " + String.valueOf(ForceZ));
//					ThreadUtil.milliSleep(500);
							          
		            while (true) {

						HumanForce = getHumanForce();

						CommandedFingerPosition = (int) (.9 * ((double) CommandedFingerPositionPrev) + .1 * (-65.0
								/ 19.0 * HumanForce + 150.0));
						CommandedFingerPositionPrev = CommandedFingerPosition;

						FingerTip.setPosReq(CommandedFingerPosition);
						// getLogger().info(String.valueOf(CommandedFingerPosition));

						if (HumanForce >= 20.0) {
							FingerTip.setPosReq(0);
							WaterBottle.detach();
							break;
						}
		            }
		            ThreadUtil.milliSleep(2000);
		            theServoRuntime.stopMotion();
					
				} else if (projectNumber == 0) { //Handover
					TcpServer t = new TcpServer();
					
					   while(true){

						  	  System.out.println("starting server");
					              t.Start();
						      System.out.println("connected");

						      t.WaitForContextRequest();
						      //estimate context, send when ready
						      double[] x = {3.1, 0.4, 0.8};
						      t.SendContext(x);
						      
//						      double[] params = t.WaitForParams();
//						      System.out.printf("First param %f\n", params[0]);
//						      //do experiment
//						      System.out.println("Experiemnt done");
						      
						      //optinal
						      //double[][] data;
						      //t.SendData(data);
						      
					             //clientSentence = t.ReadLn();
					             //System.out.println("Received: " + clientSentence);

					             //t.Write("You said:" + clientSentence + '\n');
						      break;
					   }
					
//					count +=1;
//					if (count == 100){
//						getLogger().info("Entered the handover loop 100 times");
//					}
//					
//					Frame msrPoseUpdate = new Frame();
//					msrPoseUpdate = lbr_iiwa_14_R820_1.getCurrentCartesianPosition(WaterBottle.getRootFrame());
//
//					posXmsr = msrPoseUpdate.getX();
//					posYmsr = msrPoseUpdate.getY();
//					posZmsr = msrPoseUpdate.getZ();
//					
//					RB1_x = (OptiTrackData[2] - robotBaseX)*1000;
//					RB1_y = (OptiTrackData[3] - robotBaseY)*1000;
//					RB1_z = (OptiTrackData[4] - robotBaseZ)*1000;
//
//					RB2_x = (OptiTrackData[10] - robotBaseX)*1000;
//					RB2_y = (OptiTrackData[11] - robotBaseY)*1000;
//					RB2_z = (OptiTrackData[12] - robotBaseZ)*1000;
//
//				    bodyDistance = Math.sqrt(  		Math.pow( (RB2_x), 2) 
//							+ 	Math.pow( ( RB2_y), 2) 
//							+ 	Math.pow( ( RB2_z), 2) 	)  ;
//				    
//				    handDistance = Math.sqrt( Math.pow( (RB1_x - posXmsr), 2) 
//							+ 	Math.pow( ( RB1_y - posYmsr), 2) 
//							+ 	Math.pow( ( RB1_z - posZmsr), 2) 	)  ;
//
//				    if ((bodyDistance <= 1500.0) ){//& (handDistance >= handTrackDistanceLimit)) {
//				    	followHand = true;
//				    } else {followHand = false;}
//
//				    if (RB1_z > 1000.0) {
//				    	getLogger().info("Hand raised => finish hand tracking");
//				    	theServoRuntime.stopMotion();
//				    	break;
//				    }
//
//				    if (followHand) { 
//
//				    	Frame goalFrame = new Frame();  
//				    	goalFrame = lbr_iiwa_14_R820_1.getCurrentCartesianPosition(WaterBottle.getRootFrame());
//
//				    	if ( RB1_z >= 800.0) {RB1_z = 800.0;}
//				    	if ( RB1_z <= -100.0) {RB1_z = -100.0;}
//
//				    	goalFrame.setX(RB1_x);
//				    	goalFrame.setY(RB1_y);
//				    	goalFrame.setZ(RB1_z);
//
//				    	double theta1 = Math.atan(RB1_y/RB1_x);
//				    	goalFrame.setAlphaRad(77.0*Math.PI/180.0 + theta1);
//				    	goalFrame.setBetaRad(0);
//				    	goalFrame.setGammaRad(90*Math.PI/180.0); 
//
//				    	theServoRuntime.setDestination(goalFrame, World.Current.getRootFrame());
//
//				    } 

				} else if (projectNumber == 6) {
					
					TcpServer t = new TcpServer();

					   while(true){

					  	  System.out.println("starting server");
				              t.Start();
					      System.out.println("connected");

					      t.WaitForContextRequest();
					      //estimate context, send when ready
					      double[] x = {3.1, 0.4, 0.8};
					      t.SendContext(x);
					      
//					      double[] params = t.WaitForParams();
//					      System.out.printf("First param %f\n", params[0]);
//					      //do experiment
//					      System.out.println("Experiemnt done");
					      
					      //optinal
					      //double[][] data;
					      //t.SendData(data);
					      
				             //clientSentence = t.ReadLn();
				             //System.out.println("Received: " + clientSentence);

				             //t.Write("You said:" + clientSentence + '\n');
					      break;
				          } 
					   break; } 
			  
				
			} 
		}

			catch(Exception e){
				e.printStackTrace();
			}



		}

		private int HandoverDialog() {
			
			int Choice = 99;

			while (Choice > 10) {
				Choice = getApplicationUI().displayModalDialog(
						ApplicationDialogType.QUESTION,
						"Please Select One of the Options", "Handover",
						"HandFollowing (noObject)", "ContextEstimation", "OptitrackCheck", "FixedPointHandover", "RB1 Orientation", "TCP/IP test");

				switch (Choice) {
				case 0:
					getLogger().info("Handover starting...");
					Choice = 0;
					break;
				case 1:
					getLogger().info("HandFollowing (noObject) starting...");
					getLogger().info("Move rigid body 2 within 1.5 meters and the robot will start following rigid body 1");
					Choice = 1;
					break;
				case 2:
					getLogger().info("ContextEstimation (noObject) starting...");
					Choice = 2;
					break;
				case 3:
					getLogger().info("OptitrackCheck (noObject) starting...");
					Choice = 3;
					break;
				case 4:
					getLogger().info("FixedPointHandover starting...");
					Choice = 4;
					break;
				case 5:
					getLogger().info("Orientation test");
					Choice = 5;
					break;
				case 6:
					getLogger().info("TCP/IP test");
					Choice = 6;
					break;

				default:
					break;
				}

			}

			return Choice;
		}
		
		private int TrackingDialog() {
		
			int Choice = 99;

			while (Choice > 10) {
				Choice = getApplicationUI().displayModalDialog(
						ApplicationDialogType.QUESTION,
						"Body or Hand tracking?", "Body",
						"Hand");

				switch (Choice) {
				case 0:
					getLogger().info("Body tracking starting within 1.5 meters starting...");
					Choice = 0;
					break;
				case 1:
					getLogger().info("Hand tracking (noObject) starting...");
					getLogger().info("Move rigid body 2 within 1.5 meters and the robot will start following rigid body 1");
					Choice = 1;
					break;

				default:
					break;
				}
			}

			return Choice;
		}
		
		private int DialogInitHandover() {
					
			int Choice = 99;

			while (Choice > 10) {
				Choice = getApplicationUI().displayModalDialog(
						ApplicationDialogType.QUESTION,
						"Handover Init", "Open",
						"Close","Grab Bottle and Start Handover", "Start Handover");

				switch (Choice) {
				case 0:
					getLogger().info("Opening...");
					Choice = 0;
					break;
				case 1:
					getLogger().info("Closing and attaching...");
					Choice = 1;
					break;
				case 2:
					getLogger().info("Grab and init...");
					Choice = 2;
					break;
				case 3:
					getLogger().info("Start Handover!");
					Choice = 3;
					break;

				default:
					break;
				}
			}
			return Choice;
		}
		
		private int FingerControlDialog() {
			
			int Choice = 99;

			while (Choice > 10) {
				Choice = getApplicationUI().displayModalDialog(
						ApplicationDialogType.QUESTION,
						"Finger Control test", "Open",
						"Close","Start Finger Control");

				switch (Choice) {
				case 0:
					getLogger().info("Opening...");
					Choice = 0;
					break;
				case 1:
					getLogger().info("Closing...");
					Choice = 1;
					break;
				case 2:
					getLogger().info("Starting Finger Control...");
					Choice = 2;
					break;

				default:
					break;
				}

			}

			return Choice;
		}
		
		private double getHumanForce() {

			ForceSensorData torData = lbr_iiwa_14_R820_1.getExternalForceTorque(
					ToolGripper.getFrame("/WaterBottleCP"),
					ToolGripper.getFrame("/WaterBottleCP"));

			double ForceX = torData.getForce().getX() ;
			double ForceY = torData.getForce().getY() ;
			double ForceZ = torData.getForce().getZ();
			double HumanForce = Math.sqrt(Math.pow(ForceY, 2)
					+ Math.pow(ForceZ / 6.0, 2));
			// double HumanForce = Math.sqrt(Math.pow(ForceX, 2) + Math.pow(ForceY,
			// 2) + Math.pow(ForceZ, 2));

			return HumanForce;
		}



		/**
		 * Auto-generated method stub. Do not modify the contents of this method.
		 */
		public static void main(String[] args) {
			OptitrackHandoverAndras_test app = new OptitrackHandoverAndras_test();
			app.runApplication();
		}
	}
