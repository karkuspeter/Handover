package application;

import java.util.*;
import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import static com.kuka.roboticsAPI.motionModel.BasicMotions.ptp;
import static com.kuka.roboticsAPI.motionModel.BasicMotions.spl;

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
import com.kuka.roboticsAPI.motionModel.*;
import com.kuka.roboticsAPI.motionModel.controlModeModel.CartesianImpedanceControlMode;
import com.kuka.roboticsAPI.sensorModel.ForceSensorData;
import com.kuka.roboticsAPI.uiModel.ApplicationDialogType;
import com.kuka.roboticsAPI.userInterface.ServoMotionUtilities;
import com.kuka.generated.ioAccess.*;


public class BallThrowingPeter extends RoboticsAPIApplication {
	private Controller kuka_Sunrise_Cabinet_1;
	private LBR lbr_iiwa_14_R820_1;
	private Tool ToolGripper;
	private GripperIOGroup FingerTip;
	private CartesianImpedanceControlMode mode = null ; 

	private UdpClient client = null;  //communication with OptiTrack
	private TcpServer MatlabServer;   //communication with Matlab

	Boolean _IsRunning = true;
	double[] OptiTrackData;

	
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
		getApplicationData().createFromTemplate("WaterBottle");
		
		ToolGripper.attachTo(lbr_iiwa_14_R820_1.getFlange());

		mode = new CartesianImpedanceControlMode() ; 
		mode.parametrize(CartDOF.A).setStiffness(5) ;  	// rotation about Z
		mode.parametrize(CartDOF.B).setStiffness(5);		// rotation about Y
		mode.parametrize(CartDOF.C).setStiffness(5);		// rotation about X

		mode.parametrize(CartDOF.X).setStiffness(300) ;
		mode.parametrize(CartDOF.Y).setStiffness(300) ;
		mode.parametrize(CartDOF.Z).setStiffness(300) ;
		mode.parametrize(CartDOF.ALL).setDamping(0.7) ;
		
		MatlabServer = new TcpServer(3, 30000);
	}


	public int RefreshOptitrack()
	{
		String recvStr = client.getString();
		OptiTrackData = ConvertOptitrackData(recvStr);
		return 1;
	}
	
	public double[] ConvertOptitrackData(String recvStr)
	{
		String[] optiTrackStr = recvStr.split(" ");
		double[] data= new double[optiTrackStr.length];
		for(int j=0;j<optiTrackStr.length;j++)
		{
			data[j]=Double.valueOf(optiTrackStr[j]);

		}		
		return data;
	}
	
	public double[][] ConvertRecording(List<String> buffer)
	{
		if (buffer.size() == 0)
			return new double[0][0];
		
		int dims = (ConvertOptitrackData(buffer.get(0))).length;
		double[][] tr = new double[buffer.size()][dims];
		for (int i=0; i<buffer.size(); i++)
		{
			// this may result in exception if some string contains different number of params
			tr[i] = ConvertOptitrackData(buffer.get(i));
		}
		
		return tr;
	}

	@Override 
	public void dispose(){

		_IsRunning = false ; 
		System.out.println(" closing sockets in Dispose Block"); 

		if(client != null)
		{
			client.kill();
		}
		
		if(MatlabServer != null)
		{
			MatlabServer.Close();
		}

		super.dispose();
	}


	public double[][] ThrowBall(PTP begin_frame, PTP end_frame){
		//should check first if frames are safe 
		
		getLogger().info("Moving to start position");
		lbr_iiwa_14_R820_1.move(begin_frame.setJointVelocityRel(0.5));
		
		ThreadUtil.milliSleep(500);

		// start recording tracking data
		client.StartRecording();
		
		getLogger().info("Throwing");
		lbr_iiwa_14_R820_1.move(end_frame);
		
		ThreadUtil.milliSleep(500); 
		
		// stop recording data
		List<String> buffer = client.GetRecording();
		double[][] data = ConvertRecording(buffer);
		
		getLogger().info("Moving back to start position");
		lbr_iiwa_14_R820_1.move(begin_frame.setJointVelocityRel(0.5));
		
		return data;
	}
	
	public void run() {
		try{
			MatlabServer.Start();
			getLogger().info("Connection established");
		}
		catch (IOException ex)
		{
			getLogger().error("Cannot start Tcp server for Matlab");
			return;
		}
		
		// it should all user to grab the cup first
		
		
		PTP begin_frame = ptp(
				Math.toRadians(0), 
				Math.toRadians(-57),
				Math.toRadians(0),
				Math.toRadians(-77),
				Math.toRadians(0),
				Math.toRadians(-85),
				Math.toRadians(55));

		PTP end_frame = ptp(
				Math.toRadians(0), 
				Math.toRadians(8),
				Math.toRadians(0),
				Math.toRadians(-88),
				Math.toRadians(0),
				Math.toRadians(-13),
				Math.toRadians(55));
		
		while (true)
		{
		
			// wait for params from Matlab
			double params[];
			double trajectory[][];
			try{
				getLogger().info("Wait for params");
				params = MatlabServer.WaitForParams();			
			}
			catch (IOException ex){
				getLogger().error("Cannot receive params from Matlab");
				return;			
			}
	
			getLogger().info(String.format("Params: %f; %f; %f;", params[0], params[1], params[2]));
			// ignore these and use default frames for now
			
			// execute throw
			trajectory = ThrowBall(begin_frame, end_frame);
			if (trajectory.length == 0){
				getLogger().error("Throw failed");
			}

			// send back trajectory to matlab
			try{
				MatlabServer.SendData(trajectory);		
			}
			catch (IOException ex){
				getLogger().error("Cannot send trajectory to Matlab");
				return;			
			}			
			
			// ask user
			int choice = getApplicationUI().displayModalDialog(
					ApplicationDialogType.QUESTION,
					"Do another one?", 
					"Yes", "Exit");
			if (choice != 0)
				break;

		}
		
//		ToolGripper.getFrame("/WaterBottleCP").move(ptp(getApplicationData().getFrame("/ThrowEndPos")));		
//		Spline mySpline;
//		mySpline = new Spline( spl(getFrame("/ThrowStartPos")).setCartVelocity(100),
//				spl(getFrame("/ThrowEndPos")) );
//		ToolGripper.getFrame("/WaterBottleCP").move(mySpline);


	}

	/**
	 * Auto-generated method stub. Do not modify the contents of this method.
	 */
	public static void main(String[] args) {
		BallThrowingPeter app = new BallThrowingPeter();
		app.runApplication();
	}
}
