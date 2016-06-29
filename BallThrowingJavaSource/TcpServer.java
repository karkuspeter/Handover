package application;
import java.io.*; 
import java.net.*; 

class TcpServer {    

private ServerSocket welcomeSocket;
private Socket connectionSocket;
private  BufferedReader inFromClient;
private DataOutputStream outToClient;
private int params_count = 3;
private int port = 30000;

public TcpServer(){
}

public TcpServer(int _params_count, int _port){
	params_count = _params_count;
	port = _port;
}

public Boolean Start() throws IOException{
	if (welcomeSocket == null){
	      System.out.println("welcome");
 		welcomeSocket = new ServerSocket(port);
	}
 	System.out.println("accept");
	connectionSocket = welcomeSocket.accept();

	System.out.println("accepted");
	
    inFromClient =  new BufferedReader(new InputStreamReader(connectionSocket.getInputStream()));
    outToClient = new DataOutputStream(connectionSocket.getOutputStream());

	return true;
}

public void Close(){
	try{
		if (connectionSocket != null)
			connectionSocket.close();
		if (welcomeSocket != null)
			welcomeSocket.close();
	} 
	catch (IOException ex)
	{}
}

public String ReadLn() throws IOException{
	return inFromClient.readLine();
}

public void Write(String s) throws IOException{
	outToClient.writeBytes(s);
}

/*public double[] ReadHandoverParams() throws Exception{
	String msg = inFromClient.readLine();
	double[] x = new double[1];
	x[0] = 1;
	if (msg == "Params\n")
	{
		System.out.println("ok "+msg);
		
		DataInputStream dis = new DataInputStream(connectionSocket.getInputStream());
        System.out.println(dis.readDouble()); // prints out garbage
	}
	else
	{
		System.out.println("error "+msg);
		msg = inFromClient.readLine();
		System.out.println("error "+msg);
	}
	return x;
}*/


public void SendData(double[][] data) throws IOException{
	Write("Data\n");
	Write(String.valueOf(data.length)+"\n");
	Write(String.valueOf(data[0].length)+"\n");
	for (double[] row:data){
		for (double x:row){
			Write(String.valueOf(x)+"\n");
		}
	}
	Write("End\n");
	System.out.println("Context sent");	

}

public void SendContext(double[] params) throws IOException{
	Write("Context\n");
	for (double x:params){
		Write(String.valueOf(x)+"\n");
	}
	Write("End\n");
	System.out.println("Context sent");	
}

public boolean WaitForContextRequest() throws IOException{
	String msg;

	while(true){
		System.out.println("waiting for request");
		msg = inFromClient.readLine();
		
		if(msg.equals("Context"))
		{
			System.out.println("Context request");
			return true;
		}
		else if(msg.equals("")){
		}
		else
		{
			System.out.println("error "+msg);
			return false;
		}
	}
}

public double[] WaitForParams() throws IOException{
	double[] x = new double[params_count];
	String msg;

	while(true){
		msg = inFromClient.readLine();

		if (msg.equals("Run"))
		{
			//System.out.println("getting param");
	
			if(Integer.parseInt(inFromClient.readLine()) != params_count){
				System.out.printf("wrong number of params\n");
				return new double[0];
			}
	
			for(int i=0; i<params_count; i++ ){
				x[i] = Double.parseDouble(inFromClient.readLine());
				System.out.printf("%f\n", x[i]);
			}
			
			msg = inFromClient.readLine();
			if (!msg.equals("End")){
				System.out.println("wrong command"+msg);
				return new double[0];
			} else {
				return x;
			}
		}
		else if(msg.equals("")){
		}
		else
		{
			System.out.println("error "+msg);
			return new double[0];
		}
	}
}

public static void example_main(String argv[]) throws Exception       {
	   TcpServer t = new TcpServer();

	   while(true){

	  	  System.out.println("starting server");
              t.Start();
	      System.out.println("connected");

	      t.WaitForContextRequest();
	      //estimate context, send when ready
	      double[] x = {3.1, 0.4, 0.8};
	      t.SendContext(x);
	      
	      double[] params = t.WaitForParams();
	      System.out.printf("First param %f\n", params[0]);
	      //do experiment
	      System.out.println("Experiemnt done");
	      
	      //optinal
	      //double[][] data;
	      //t.SendData(data);
	      
             //clientSentence = t.ReadLn();
             //System.out.println("Received: " + clientSentence);

             //t.Write("You said:" + clientSentence + '\n');
          }
       }
 }