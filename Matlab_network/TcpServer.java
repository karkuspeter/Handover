import java.io.*; 
import java.net.*; 
class TcpServer {    

private ServerSocket welcomeSocket;
private Socket connectionSocket;
private  BufferedReader inFromClient;
private DataOutputStream outToClient;
private int params_count = 3;

public TcpServer(){
}

public Boolean Start() throws Exception{
	if (welcomeSocket == null){
	      System.out.println("welcome");
 		welcomeSocket = new ServerSocket(6788);
	}
 	System.out.println("accept");
	connectionSocket = welcomeSocket.accept();

       inFromClient =  new BufferedReader(new InputStreamReader(connectionSocket.getInputStream()));
       outToClient = new DataOutputStream(connectionSocket.getOutputStream());

	return true;
}

public String ReadLn() throws Exception{
	return inFromClient.readLine();
}

public void Write(String s) throws Exception{
	outToClient.writeBytes(s);
}

public double[] ReadHandoverParams() throws Exception{
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
}

public double[] WaitForCommand() throws Exception{
	DataInputStream dis = new DataInputStream(connectionSocket.getInputStream());
	DataOutputStream dos = new DataOutputStream(connectionSocket.getOutputStream());

	double[] x = new double[params_count];
	String msg;

	while(true){
	msg = inFromClient.readLine();

	if (msg.equals("Run"))
	{
		System.out.println("getting param");

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
	else if(msg.equals("Context"))
	{
		System.out.println("Context request");
		Write("Context\n");
		Write(String.valueOf(3.14)+"\n");
		Write(String.valueOf(3.16)+"\n");
		Write(String.valueOf(3.18)+"\n");
		Write("End\n");
		System.out.println("Context sent");	
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

public static void main(String argv[]) throws Exception       {
          String clientSentence;
          String capitalizedSentence;
	   TcpServer t = new TcpServer();

	   while(true){

	  	  System.out.println("starting server");
              t.Start();
	      System.out.println("connected");

	 	//t.Write("Hello gyurcsany\n");
		t.WaitForCommand();
		System.out.println("Experiemnt done");
             //clientSentence = t.ReadLn();
             //System.out.println("Received: " + clientSentence);

             //t.Write("You said:" + clientSentence + '\n');
          }
       }
 }